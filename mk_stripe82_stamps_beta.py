from astropy.io import fits
from astropy.cosmology.realizations import Planck15 as cosmo
import numpy as np
import matplotlib.pyplot as plt
import os
import gzip
import shutil
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

    ang = np.radians(pa)# + 90.)
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
        im[i,:] = np.sqrt((ytemp*ratio)**2 + xtemp**2)

    return im

def compress_gzip(file, compresslevel=6):
    with open(file, 'rb') as f_in:
        with gzip.open(f'{file}.gz', mode='wb', compresslevel=compresslevel) as f_comp:
            shutil.copyfileobj(f_in, f_comp)
    os.remove(file)


def get_hdu(fname):
    hdu_side = fits.open(indir+"f"+fname+".rec.fits.gz")
    ra = hdu_side[0].header["CRVAL1"]
    ra_pix = hdu_side[0].header["CRPIX1"]
    dec = hdu_side[0].header["CRVAL2"]
    dec_pix = hdu_side[0].header["CRPIX2"]
    data = hdu_side[0].data
    return ra, ra_pix, dec, dec_pix, data


with open("config.yaml", "r") as f: 
    config = yaml.load_all(f)

indir = config["obs"]["stripe82_frames_path"]
outdir = config["obs"]["out_path"]
plotdir = outdir + "jpgs/"

hdu_cat = fits.open(config["obs"]["catalogue"])

fname = hdu_cat[1].data["fname"]
gal_id = hdu_cat[1].data["JID"]
ra_gal = hdu_cat[1].data["RAJ2000"]#["RA"]
dec_gal = hdu_cat[1].data["DEJ2000"]#["DEC"]
Redge_deg = hdu_cat[1].data["Redge"]*cosmo.arcsec_per_kpc_proper(hdu_cat[1].data["z"]).value/3600.
PA = hdu_cat[1].data["PA"]
q_ellip = hdu_cat[1].data["q"]

mid_size_pix = np.floor(config["obs"]["cutout_size"]/2)


sufix = "_stamps"+str(mid_size_pix*2)+".fits"
missing_or_corrupted = []

if not os.path.exists(outdir+"jpgs"):
    sys("mkdir "+outdir+"jpgs")

for ii in range(len(fname)):
    print("Ngal:", ii)
    padded = False
    try:
        print("###########################################")
        if os.path.exists(outdir+gal_id[ii]+sufix+".gz"):
            continue

        stamps = np.zeros((mid_size_pix*2, mid_size_pix*2, 4), dtype=np.single)
        masks = np.zeros((mid_size_pix*2, mid_size_pix*2), dtype=np.single)

        ra_plus_name = fname[ii]
        dec_plus_name = fname[ii]

        #plot = True        

        if os.path.exists(indir+"f"+fname[ii]+"_u.rec.fits.gz"):
            hdu = fits.open(indir+"f"+fname[ii]+"_u.rec.fits.gz")
            if np.single(fname[ii][:3])>100:
                ra_cent = (np.single(fname[ii][:3])-101.) * 0.5+0.25#hdu[0].header["CRVAL1"]
            else:
                ra_20_41_deg = (20. + 41./60.) * 180./12.
                ra_cent = ra_20_41_deg + (np.single(fname[ii][:3])-1.) * 0.5
            dec_cent = (np.single(fname[ii][3])) * 0.5 -1.5#hdu[0].header["CRVAL2"]
            x_cent = hdu[0].header["CRPIX1"]
            y_cent = hdu[0].header["CRPIX2"]
            ra_pix_scale = hdu[0].header["CDELT1"]
            dec_pix_scale = hdu[0].header["CDELT2"]

            if ii==0:
                plt.hist(Redge_deg/dec_pix_scale, bins=50);
                plt.ylabel("N galaxies")
                plt.xlabel(r"R$_{edge}$ [pix]")
                plt.xlim(0,470)
                plt.axvline(300, c="r", ls="--")
                plt.axvline(200, c="r", ls="--")
                plt.axvline(400, c="r", ls="--")
                plt.text(70, 50, "total galaxies: "+str(Redge_deg.size))
                plt.text(205, 25, "400 pix")
                plt.text(305, 25, "600 pix")
                plt.text(405, 25, "800 pix")
                plt.text(205, 20, str(np.array([Redge_deg/dec_pix_scale<200]).sum())+" gal")
                plt.text(305, 20, str(np.array([Redge_deg/dec_pix_scale<300]).sum())+" gal")
                plt.text(405, 20, str(np.array([Redge_deg/dec_pix_scale<400]).sum())+" gal")
                plt.savefig(plotdir+"hist_regde_pix.png", bbox_inches='tight')
                plt.close()



            hdu_shape = hdu[0].data.shape
            half_data = np.int32((hdu_shape[0]-6)/2)

            x_gal = np.int32(x_cent + (ra_gal[ii]-ra_cent)/ra_pix_scale)
            y_gal = np.int32(y_cent + (dec_gal[ii]-dec_cent)/dec_pix_scale)


            x_gal_final = -999
            y_gal_final = -999
            fits_in_frame = True

            fits_in_frame_ra = (x_gal-mid_size_pix)<0 or (x_gal+mid_size_pix)>=hdu_shape[0]
            
            fits_in_frame_dec = ((y_gal-mid_size_pix)<0 or (y_gal+mid_size_pix)>=hdu_shape[1])
            

        else:
            missing_or_corrupted.append(fname[ii])
            print("stamp doesn't fit, try downloading more frames with download_update.sh")
            continue

        if fits_in_frame_ra and fits_in_frame_dec:

            ra_plus_name = str(np.int32(fname[ii][:3])- np.sign(x_gal-mid_size_pix)).zfill(3) + fname[ii][-1]
            dec_plus_name = str(np.int32(fname[ii][:4])+np.sign(y_gal-mid_size_pix)).zfill(4)
            radec_plus_name = (str(np.int32(fname[ii][:3])- np.sign(x_gal-mid_size_pix))+\
                            str(np.int32(fname[ii][3])+np.sign(y_gal-mid_size_pix))).zfill(4)
            

            if (ra_plus_name[:-1]=="100") or (ra_plus_name[:-1]=="220"):
                ra_plus_name = "zero"
                radec_plus_name = "zero"
                padded = True
                #print("stamp doesn't fit in observations")
                #continue

            if dec_plus_name[-1]=="0" or dec_plus_name[-1]=="6":
                dec_plus_name = "zero"
                radec_plus_name = "zero"
                padded = True
                #print("stamp doesn't fit in observations")
                #continue

            needed_frames = True

            if not os.path.exists(indir+"f"+ra_plus_name+"_g.rec.fits.gz"):
                missing_or_corrupted.append(ra_plus_name)
                print("stamp doesn't exist, try downloading more frames with download_update.sh")
                needed_frames = False
                pass

            if not os.path.exists(indir+"f"+dec_plus_name+"_g.rec.fits.gz"):
                missing_or_corrupted.append(dec_plus_name)
                print("stamp doesn't exist, try downloading more frames with download_update.sh")
                needed_frames = False
                pass

            if not os.path.exists(indir+"f"+radec_plus_name+"_g.rec.fits.gz"):
                missing_or_corrupted.append(radec_plus_name)
                print("stamp doesn't exist, try downloading more frames with download_update.sh")
                needed_frames = False
                pass

            if not needed_frames:
                continue


            # reorganize fits such that 
            # file1 is at the top left corner
            # file2 is at the top right corner
            # file3 is at the bottom left
            # file4 is at the bottom right
            

            if (x_gal-mid_size_pix)<0:
                #main image is on the right side
                if (y_gal-mid_size_pix)<0:
                    #main image is on the bottom side
                    fname1 = radec_plus_name
                    fname2 = dec_plus_name
                    fname3 = ra_plus_name
                    fname4 = fname[ii]
                    central_pix = [half_data*3, half_data*3]

                else:
                    #main image is on the top side
                    fname1 = ra_plus_name
                    fname2 = fname[ii]  
                    fname3 = radec_plus_name
                    fname4 = dec_plus_name
                    central_pix = [half_data*3, half_data]  

            else:    
                #main image is on the left side  
                if (y_gal-mid_size_pix)<0:
                    #main image is on the bottom side
                    fname1 = dec_plus_name  
                    fname2 = radec_plus_name
                    fname3 = fname[ii]
                    fname4 = ra_plus_name
                    central_pix = [half_data, half_data*3] 
                    
                else:
                    #main image is on the top side
                    fname1 = fname[ii]
                    fname2 = ra_plus_name
                    fname3 = dec_plus_name
                    fname4 = radec_plus_name
                    central_pix = [half_data, half_data] 
                            

            megastamp = np.zeros((4*half_data, 4*half_data, 5), dtype=np.single)
            for jj, band in enumerate(["u", "g", "r", "i", "z"]):
                megastamp[:2*half_data, :2*half_data, jj] = \
                    fits.getdata(indir+"f"+fname1+"_"+band+".rec.fits.gz")[4:-3, 4:-3]
                megastamp[:2*half_data, 2*half_data:, jj] = \
                    fits.getdata(indir+"f"+fname2+"_"+band+".rec.fits.gz")[4:-3, 4:-3]
                megastamp[2*half_data:, :2*half_data, jj] = \
                    fits.getdata(indir+"f"+fname3+"_"+band+".rec.fits.gz")[4:-3, 4:-3]
                megastamp[2*half_data:, 2*half_data:, jj] = \
                    fits.getdata(indir+"f"+fname4+"_"+band+".rec.fits.gz")[4:-3, 4:-3]
            
        elif fits_in_frame_ra:

            ra_plus_name = str(np.int32(fname[ii][:3])- np.sign(x_gal-mid_size_pix)).zfill(3) + fname[ii][-1]

            if (ra_plus_name[:-1]=="100") or (ra_plus_name[:-1]=="220"):
                ra_plus_name = "zero"
                padded = True


            if not os.path.exists(indir+"f"+ra_plus_name+"_g.rec.fits.gz"):
                missing_or_corrupted.append(ra_plus_name)
                print("stamp doesn't exist, try downloading more frames with download_update.sh")
                continue


            if (x_gal-mid_size_pix)<0:
                #main image is on the right side
                fname1 = ra_plus_name
                fname2 = fname[ii]
                central_pix = [half_data*3, half_data]
            else:
                fname1 = fname[ii]
                fname2 = ra_plus_name
                central_pix = [half_data,half_data]

            megastamp = np.zeros((2*half_data, 4*half_data, 5), dtype=np.single)
            for jj, band in enumerate(["u", "g", "r", "i", "z"]):
                megastamp[:, :2*half_data, jj] = \
                    fits.getdata(indir+"f"+fname1+"_"+band+".rec.fits.gz")[4:-3, 4:-3]
                megastamp[:, 2*half_data:, jj] = \
                    fits.getdata(indir+"f"+fname2+"_"+band+".rec.fits.gz")[4:-3, 4:-3]
        
        
        elif fits_in_frame_dec:

            dec_plus_name = str(np.int32(fname[ii][:4])+np.sign(y_gal-mid_size_pix)).zfill(4)

            if dec_plus_name[-1]=="0" or dec_plus_name[-1]=="6":
                dec_plus_name = "zero"
                padded = True
                #print("stamp doesn't fit in observations")
                #continue

            if not os.path.exists(indir+"f"+dec_plus_name+"_g.rec.fits.gz"):
                missing_or_corrupted.append(dec_plus_name)
                print("stamp doesn't exist, try downloading more frames with download_update.sh")
                continue

            if (y_gal-mid_size_pix)<0:
                #main image is on the bottom
                fname1 = dec_plus_name
                fname2 = fname[ii]
                central_pix = [half_data, half_data*3]
            else:
                fname1 = fname[ii]
                fname2 = dec_plus_name
                central_pix = [half_data, half_data]

            megastamp = np.zeros((4*half_data, 2*half_data, 5), dtype=np.single)
            for jj, band in enumerate(["u", "g", "r", "i", "z"]):
                megastamp[:2*half_data, :, jj] = \
                    fits.getdata(indir+"f"+fname1+"_"+band+".rec.fits.gz")[4:-3, 4:-3]
                megastamp[2*half_data:, :, jj] = \
                    fits.getdata(indir+"f"+fname2+"_"+band+".rec.fits.gz")[4:-3, 4:-3]

        else:

            megastamp = np.zeros((hdu_shape[0], hdu_shape[1], 5), dtype=np.single)
            for jj, band in enumerate(["u", "g", "r", "i", "z"]):
                megastamp[:, :, jj] = \
                    fits.getdata(indir+"f"+fname[ii]+"_"+band+".rec.fits.gz")
            central_pix = [x_cent, y_cent]

        x_gal = np.int32(central_pix[0] + (ra_gal[ii]-ra_cent)/ra_pix_scale)
        y_gal = np.int32(central_pix[1] + (dec_gal[ii]-dec_cent)/dec_pix_scale)

        #print(central_pix)
        #print(x_gal, y_gal)
        #print((ra_gal[ii]-ra_cent)/ra_pix_scale, (dec_gal[ii]-dec_cent)/dec_pix_scale)
        chamba_stamp = megastamp[y_gal-mid_size_pix:y_gal+mid_size_pix, \
                                 x_gal-mid_size_pix:x_gal+mid_size_pix, :]
        
        #if plot:
        #    plt.imshow(np.log10(megastamp[:,:,0]))
        #    plt.plot(central_pix[0], central_pix[1], "r", ms=10, ls='', marker='o', alpha=0.5)
        #    plt.plot(x_gal, y_gal, "r", ms=10, ls='', marker='o', alpha=0.2)
        #    plt.show()

        #    plt.imshow(np.log10(chamba_stamp[:,:,0]))
        #    plt.show()

        masks = dist_ellipse([2*mid_size_pix, 2*mid_size_pix], mid_size_pix+0.5, 
                                 mid_size_pix+0.5, (1./q_ellip[ii]), pa=PA[ii])
        fake_color = np.zeros((chamba_stamp.shape[0], chamba_stamp.shape[1], 3))
        
        for jj in range(1,4):
            fake_color[:,:,jj-1] = np.log10(chamba_stamp[:,:,jj]-chamba_stamp[:,:,jj].min())/np.max(np.log10(chamba_stamp[:,:,jj]-chamba_stamp[:,:,jj].min()))
        plt.imshow(fake_color)
        plt.imshow(np.where(masks<np.abs(Redge_deg[ii]/ra_pix_scale), 1,0), alpha=0.1, cmap='bone')
        plt.tick_params(left = False, right = False , labelleft = False , 
            labelbottom = False, bottom = False) 
        plt.tight_layout()
        plt.savefig(outdir+"jpgs/"+gal_id[ii]+".jpg")
        plt.close()
    
        primary = fits.PrimaryHDU(chamba_stamp)
        mask_ext = fits.ImageHDU(masks)

        h = primary.header
        h["CRVAL1"] = ra_gal[ii]
        h["CRPIX1"] = mid_size_pix
        h["CDELT1"] = ra_pix_scale
        h["CRVAL2"] = dec_gal[ii]
        h["CRPIX2"] = mid_size_pix
        h["CDELT2"] = dec_pix_scale
        h["PADDED"] = padded

        h = mask_ext.header
        h["EXTNAME"] = "MASK"
        h["REDGE"] = np.abs(Redge_deg[ii]/ra_pix_scale)
        h["REDGE_U"] = "[pixels]"

        hlist = fits.HDUList([primary, mask_ext])
        hlist.update_extend()
        out_fit = outdir+gal_id[ii]+sufix
        hlist.writeto(out_fit, overwrite=1)
        compress_gzip(out_fit)


    except KeyError:
        print("CRVAL missing")


missing_or_corrupted = np.intersect1d(missing_or_corrupted, missing_or_corrupted)
print("Missing or corrupted:", len(missing_or_corrupted))

f = open(config["obs"]["download_script"], "w")
f.write("cd "+ indir +"\n")
for ii in missing_or_corrupted:
    f.write("curl -O ftp://stripero:s0l0le0@ftp.iac.es/coadds/f"+ii+"_u.rec.fits\n")
    f.write("gzip f"+ii+"_u.rec.fits\n")
    f.write("curl -O ftp://stripero:s0l0le0@ftp.iac.es/coadds/f"+ii+"_g.rec.fits\n")
    f.write("gzip f"+ii+"_g.rec.fits\n")
    f.write("curl -O ftp://stripero:s0l0le0@ftp.iac.es/coadds/f"+ii+"_r.rec.fits\n")
    f.write("gzip f"+ii+"_r.rec.fits\n")
    f.write("curl -O ftp://stripero:s0l0le0@ftp.iac.es/coadds/f"+ii+"_i.rec.fits\n")
    f.write("gzip f"+ii+"_i.rec.fits\n")
    f.write("curl -O ftp://stripero:s0l0le0@ftp.iac.es/coadds/f"+ii+"_z.rec.fits\n")
    f.write("gzip f"+ii+"_z.rec.fits\n")

