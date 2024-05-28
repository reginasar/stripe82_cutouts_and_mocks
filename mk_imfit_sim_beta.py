import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.cosmology.realizations import Planck15 as cosmo
import matplotlib.pyplot as plt
import os
import gzip
import shutil
import yaml

with open("config.yaml", "r") as f: 
    config = yaml.load_all(f)


def compress_gzip(file, compresslevel=6):
    with open(file, 'rb') as f_in:
        with gzip.open(f'{file}.gz', mode='wb', compresslevel=compresslevel) as f_comp:
            shutil.copyfileobj(f_in, f_comp)
    os.remove(file)

def mk_config_file(info_dict, config_name):
    with open(config_name, "w") as f:
        for params in info_dict.items():
            f.write(params[0]+" "+str(params[1])+"\n")


bands = ["u", "g", "r", "i", "z"]

psf_file = ["--psf="+config["sim"]["psf_file"].replace("?", band)+\
            " " for band in bands]
            
pix_scale = 0.396 * u.arcsec
img_size = np.floor(config["sim"]["cutout_size"])
centre_pix = img_size/2 + 0.5

hdu_cat = fits.open(config["sim"]["catalogue"])

n_gal = hdu_cat[1].data["PA_disc"].size

rng = np.random.default_rng(12345)

z_gal = hdu_cat[1].data["z_gal"]
pa_disc = hdu_cat[1].data["PA_disc"]
incl_disc = hdu_cat[1].data["incl_disc"]
r_edge = hdu_cat[1].data["r_edge"]* u.kpc * cosmo.arcsec_per_kpc_proper(z_gal) / pix_scale
h1_ = hdu_cat[1].data["h1_"]* u.kpc * cosmo.arcsec_per_kpc_proper(z_gal) / pix_scale
h2_ = hdu_cat[1].data["h2_"]* u.kpc * cosmo.arcsec_per_kpc_proper(z_gal) / pix_scale
alpha = hdu_cat[1].data["alpha"]/ u.kpc / cosmo.arcsec_per_kpc_proper(z_gal) * pix_scale
z0 = hdu_cat[1].data["z0"]* u.kpc * cosmo.arcsec_per_kpc_proper(z_gal) / pix_scale
n_disc = hdu_cat[1].data["n_disc"]
PA_bulge = hdu_cat[1].data["PA_bulge"]
ell_bulge = hdu_cat[1].data["ell_bulge"]
Re_bulge = hdu_cat[1].data["Re_bulge"]* u.kpc * cosmo.arcsec_per_kpc_proper(z_gal) / pix_scale
n_bulge = hdu_cat[1].data["n_bulge"]
I_0 = 1.
J_0 = I_0/(2*z0)


outdir = config["sim"]["config_files_path"]
outdir_fits = config["sim"]["out_path"]
f = open(config["sim"]["gen_file"], "w")
f.write("#!/bin/bash\n")
f.write("\n")
f.write("MKIMG="+config["sim"]["make_image"])
f.write("\n")
f.write("cd "+outdir+"\n")
f.write("\n")


position = {"X0": centre_pix,
            "Y0": centre_pix,
            }



for ii in range(1000):

    ########################################################################
    ######## BULGES ########################################################
    ########################################################################

    bulge = {"FUNCTION": "Sersic # LABEL bulge",
             "PA": PA_bulge[ii],
             "ell": ell_bulge[ii],
             "n": n_bulge[ii],
             "I_e": J_0[ii]/10.,
             "r_e": Re_bulge[ii],
             }
    final_dict = position.copy()
    final_dict.update(bulge)
    name = "bulge_"+str(ii)

    mk_config_file(final_dict, outdir+name+".dat")

    for jj,band in enumerate(bands):
        f.write("$MKIMG --ncols="+str(img_size)+" --nrows="+str(img_size)+\
                " "+psf_file[jj]+"-o="+outdir_fits+name+"_"+band+".fits "+\
                name+".dat\n")
        f.write("gzip -f "+outdir_fits+name+"_"+band+".fits\n")


    ########################################################################
    ######## TRUNCATED DISC ################################################
    ########################################################################


    trunc_disc = {"FUNCTION": "BrokenExponentialDisk3D # LABEL truncation_disc",
             "PA": pa_disc[ii],    # Position angle ,
             "inc": incl_disc[ii], # cos^-1 (b/a),
             "J_0": J_0[ii],  # Central luminosity density (on area),
             "h1": h1_[ii],     # Scale length of the disc (pix),
             "h2": h2_[ii],   # 0,80  Scale length halo (pix),
             "r_break": r_edge[ii],   # Break radius,
             "alpha": alpha[ii],  # 0.1, 1The smaller the smother the break,
             "n": n_disc[ii],    # the higher the more similar to an exp.,
             "z_0": z0[ii],    # Vertical scale height,
             }

    final_dict = position.copy()
    final_dict.update(trunc_disc)
    name = "trunc_"+str(ii)

    mk_config_file(final_dict, outdir+name+".dat")

    for jj,band in enumerate(bands):
        f.write("$MKIMG --ncols="+str(img_size)+" --nrows="+str(img_size)+\
            " "+psf_file[jj]+"-o="+outdir_fits+name+"_"+band+".fits "+name+\
            ".dat\n")
        f.write("gzip -f "+outdir_fits+name+"_"+band+".fits\n")


f.close()



