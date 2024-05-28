import numpy as np
from astropy.io import fits
import os
import gzip
import shutil
import astropy.units as u
from astropy.cosmology.realizations import Planck15 as cosmo
import glob
import matplotlib.pyplot as plt
import yaml
import sys


def dist_ellipse(n, xc, yc, ratio, pa=0): 
    """
    original implementation (like DIST_ELLIPSE IDL function)
    
    N = either  a scalar specifying the size of the N x N square output
              array, or a 2 element vector specifying the size of the
               M x N rectangular output array.
       XC,YC - Scalars giving the position of the ellipse center.   This does
               not necessarily have to be within the image
       RATIO - Scalar giving the ratio of the major to minor axis.   This
               should be greater than 1 for position angle to have its
               standard meaning.
    OPTIONAL INPUTS:
      POS_ANG - Position angle of the major axis in degrees, measured counter-clockwise
               from the Y axis.  For an image in standard orientation
               (North up, East left) this is the astronomical position angle.
               Default is 0 degrees.
    OUTPUT:
       IM - REAL*4 elliptical mask array, of size M x N.  THe value of each
               pixel is equal to the semi-major axis of the ellipse of center
                XC,YC, axial ratio RATIO, and position angle POS_ANG, which
               passes through the pixel.
    """

    ang = np.radians(pa + 90.)
    cosang = np.cos(ang)
    sinang = np.sin(ang)
    nx = n[1]
    ny = n[0]
    x = np.arange(-xc,nx-xc)
    y = np.arange(-yc,ny-yc)

    im = np.empty(n)
    xcosang = x*cosang
    xsinang = x*sinang

    for i in range(0, ny):

        xtemp = xcosang + y[i]*sinang
        ytemp = -xsinang + y[i]*cosang
        im[i,:] = np.sqrt((ytemp*ratio)**2 + (xtemp)**2)

    return im

def compress_gzip(file, compresslevel=6):
    with open(file, 'rb') as f_in:
        with gzip.open(f'{file}.gz', mode='wb', compresslevel=compresslevel) as f_comp:
            shutil.copyfileobj(f_in, f_comp)
    os.remove(file)



with open("config.yaml", "r") as f: 
    config = yaml.load_all(f)

indir = config["sim"]["simulated_files_path"]
outdir = config["sim"]["out_path"]

rng = np.random.default_rng(12345)

img_size = config["sim"]["cutout_size"]
pix_scale = 0.396 * u.arcsec
margin_bg = 10

hdu_cat = fits.open(config["sim"]["catalogue"])

sed_ind = hdu_cat[1].data["sed_ind"]
sed_factor = hdu_cat[1].data["sed_factor"]
bulge_fraction = hdu_cat[1].data["fraction_b"]
z_gal = hdu_cat[1].data["z_gal"]
PA = hdu_cat[1].data["pa_disc"]
incl = hdu_cat[1].data["incl_disc"] * np.pi / 180.
z0_pix = hdu_cat[1].data["z0"] * u.kpc * cosmo.arcsec_per_kpc_proper(z_gal) / pix_scale
major_axis_pix = hdu_cat[1].data["r_edge"] * u.kpc * cosmo.arcsec_per_kpc_proper(z_gal) / pix_scale
minor_axis_pix = np.array([np.max((np.abs(major_axis_pix[ii].value * np.cos(incl[ii])),
                                  np.abs(z0_pix[ii].value * np.sin(incl[ii])))) for ii in range(len(z_gal))])
q_ellip =  major_axis_pix / minor_axis_pix 


seds_20 = np.load(config["obs"]["seds_path"]+"seds_path.npy")

if not os.path.exists(outdir+"jpgs"):
    sys("mkdir "+outdir+"jpgs")

for ii in range(1000):
    ##################################################
    ##### RESCALE AND STACK FITS #####################
    ##################################################
    
    mock_gal = np.zeros((img_size, img_size, 5), dtype=np.single)
    mock_bulge = np.zeros_like(mock_gal)
    mock_disc = np.zeros_like(mock_gal)

    target_sed = seds_20[sed_ind[ii]] * sed_factor[ii]

    fraction_b = bulge_fraction[ii]


    if not os.path.exists(indir+"bulge_"+str(ii)+"_g.fits.gz"):
        continue
    
    for jj, band in enumerate(["u","g","r","i","z"]):

        bulge_out, b_head = fits.getdata(indir+"bulge_"+str(ii)+"_"+band+".fits.gz", header=True)
        mock_bulge[:,:,jj] = bulge_out

        trunc_out, t_head = fits.getdata(indir+"trunc_"+str(ii)+"_"+band+".fits.gz", header=True)
        mock_disc[:,:,jj] = trunc_out

        flux_b_original = np.nansum(bulge_out)
        flux_d_original = np.nansum(trunc_out)

        flux_b_update = bulge_out * fraction_b /flux_b_original
        flux_d_update = trunc_out * (1.-fraction_b) /flux_d_original

        flux_img = flux_b_update + flux_d_update
        xmin_ = - np.floor_divide(config["sim"]["cutout_size"], 2) + 0.5
        xmax_ = - xmin_ + 0.1
        X,Y = np.mgrid[xmin_:xmax_:1, xmin_:xmax_:1]

        apper_mask = np.where((X**2 + Y**2 < 400), 1, 0)
        flux_img_20 = np.nansum(flux_img*apper_mask)        
        
        flux_img_final = flux_img * target_sed[jj] / flux_img_20 
        mock_gal[:,:,jj] += flux_img_final


    ##################################################
    ##### ADD SDSS BACKGROUND ########################
    ##################################################

    bg_frame_name = rng.choice(glob.glob(config["sim"]["background_frames_path"]+"*_u.rec.fits.gz"))
    hdu, header = fits.getdata(bg_frame_name, header=True)

    frame_min, frame_max = margin_bg, hdu.shape[0]- margin_bg - img_size
    xmin, ymin = tuple(rng.choice(np.arange(frame_min, frame_max), 2))
    xmax, ymax = xmin+img_size, ymin+img_size
    mock_gal[:,:,0] += hdu[xmin:xmax, ymin:ymax]
    

    for jj, band in enumerate(["g", "r", "i", "z"]):
        print(bg_frame_name[:-13] + band + bg_frame_name[-12:])
        bg_frame_name = bg_frame_name[:-13] + band + bg_frame_name[-12:]
        hdu = fits.getdata(bg_frame_name)
        mock_gal[:,:,jj+1] += hdu[xmin:xmax, ymin:ymax]


    ####################################################
    ########### MAKE MASK ##############################
    ####################################################

    masks = dist_ellipse([img_size, img_size], img_size/2.,
                                 img_size/2., (q_ellip[ii]), pa=PA[ii])
    masks = np.array(masks, dtype=np.single)

    ###################################################
    ########## MAKE FITS OUTPUT FILE ##################
    ###################################################

    h0 = fits.PrimaryHDU(np.transpose(mock_gal, (2,0,1)))
    h = h0.header
    h["FRAC_BUL"] = fraction_b

    try:
        h["RA_SDSS"] = header["CRVAL1"] + (xmin+img_size/2) * header["CDELT1"]
        h["DEC_SDSS"] = header["CRVAL2"] + (ymin+img_size/2) * header["CDELT2"]
    except:
        print("no info in header")
    h["BG_FNAME"] = bg_frame_name[-18:-14]
    h["BG_XMIN"] = xmin
    h["BG_XMAX"] = xmax
    h["BG_YMIN"] = ymin
    h["BG_YMAX"] = ymax
    
    h1 = fits.ImageHDU(np.transpose(mock_bulge, (2,0,1)), b_head, "BULGE")
    h2 = fits.ImageHDU(np.transpose(mock_disc, (2,0,1)), t_head, "TRUNCATEDDISC")
    h3 = fits.ImageHDU(target_sed, name="SED_target")
    h4 = fits.ImageHDU(seds_20[sed_ind[ii]], name="SED_original")
    h5 = fits.ImageHDU(masks, name="MASK")
    hh = h5.header
    hh["REDGE"] = np.abs(major_axis_pix[ii].value)
    hh["REDGE_U"] = "[pixels]"

    
    hlist = fits.HDUList([h0, h1, h2, h3, h4, h5])
    hlist.update_extend()
    out_fit = outdir+"SDSSlike_"+str(ii)+".fits"
    hlist.writeto(out_fit, overwrite=1)
    compress_gzip(out_fit)

    ############################################
    ########## MAKE JPG ########################
    ############################################

    fake_color = np.zeros((mock_gal.shape[0], mock_gal.shape[1], 3))
        
    for jj in range(1,4):
        fake_color[:,:,jj-1] = np.log10(mock_gal[:,:,jj]-mock_gal[:,:,jj].min())/np.max(np.log10(mock_gal[:,:,jj]-mock_gal[:,:,jj].min()))
    
    plt.imshow(fake_color)
    plt.imshow(np.where(masks<major_axis_pix[ii].value, 1,0), alpha=0.1, cmap='bone')
    plt.text(10, 580, "gal_id "+str(ii), c="white", size=10)
    plt.text(10, 560, "BG-frame "+bg_frame_name[-18:-14], c="white", size=10)
    plt.text(10, 540, "PA "+f'{PA[ii]:.1f}', c="white", size=10)
    plt.text(10, 520, "incl "+f'{hdu_cat[1].data["incl_disc"][ii]:.1f}', c="white", size=10)
    plt.text(10, 500, "bulge_fr "+f'{fraction_b:.2f}', c="white", size=10)
    plt.text(10, 480, "major semiax "+f'{major_axis_pix[ii].value:.1f}', c="white", size=10)
    plt.tick_params(left = False, right = False , labelleft = False , 
            labelbottom = False, bottom = False) 
    plt.tight_layout()
    plt.savefig(outdir+"jpgs/mock_"+str(ii)+".jpg", dpi=300)
    plt.close()

