import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import ast
import yaml

with open("config.yaml", "r") as f: 
    config = yaml.load_all(f)

stripe_path = config["obs"]["stripe82_frames_path"]
hdu = fits.open(config["obs"]["catalogue"])

xmin = - np.floor_divide(config["obs"]["cutout_size"],2) + 0.5
xmax = - xmin + 0.1
X,Y = np.mgrid[xmin:xmax:1, xmin:xmax:1]
apper_mask = np.where((X**2 + Y**2 < 400), 1, 0) # circular apperture of 20 pix radius mask
apper_mask_bands = np.transpose(np.tile(apper_mask, (5,1,1)), (1,2,0))

fsed_20 = np.zeros((hdu[1].data["JID"].shape[0], 5))

for ii, gal_id in enumerate(hdu[1].data["JID"]):

    cutout_name = config["obs"]["out_path"]+gal_id+"_stamps"+\
                      str(config["obs"]["cutout_size"])+".fits.gz"

    if os.path.exists(cutout_name):

        data = fits.getdata(cutout_name)
        fsed_20[ii] = np.sum(data*apper_mask_bands, axis=(0,1))#\ assume background subtracted
            # - (np.median(data, axis=(0,1))*apper_mask_bands.sum(axis=(0,1))) #subtract background median

np.save(config["obs"]["seds_path"]+"seds_obs", fsed_20)
    


