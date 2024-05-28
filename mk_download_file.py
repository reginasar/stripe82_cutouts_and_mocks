from astropy.io import fits
import numpy as np
import matplotlib as plt
from astropy.table import Table
import yaml

with open("config.yaml", "r") as f: 
    config = yaml.load_all(f)


hdu = fits.open(config["obs"]["catalogue"])

ra = hdu[1].data["RAJ2000"] #* 180./12. 
dec = hdu[1].data["DEJ2000"]

ra_3_57_deg = (3. + 57./60.) * 180./12. - 0.25
ra_20_41_deg = (20. + 41./60.) * 180./12.
ra_23_59_deg = (23. + 59./60.) * 180./12. 

ra_bins1 = np.linspace(0., ra_3_57_deg, 119) 
ra_bins2 = np.linspace(ra_20_41_deg+0.25, ra_23_59_deg+0.25, 100)
dec_bins = np.linspace(-0.75, 1.25, 5)

dec_ind = np.digitize(dec, dec_bins) + 1
ra_ind = np.zeros_like(dec_ind)
for ii in range(len(ra)):
    if ra[ii] < (4.*180./12.):
        ra_ind[ii] = np.digitize(ra[ii], ra_bins1) + 100
    else:
        ra_ind[ii] = np.digitize(ra[ii], ra_bins2) + 1

num_str = [str(ra_ind[ii]).zfill(3)+str(dec_ind[ii]) for ii in range(len(ra))]
limited_num = np.intersect1d(num_str, num_str)

#if 'fname' not in hdu[1].data.fields(1):
new_cols = fits.ColDefs([fits.Column(name="fname", format="4A",
                                      array=num_str)])
orig_cols = hdu[1].data.columns
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
hdu.writeto(config["obs"]["catalogue"], overwrite=True)


with open(config["obs"]["download_script"], "w") as f:
    f.write("cd "+config["obs"]["stripe82_frames_path"]+"\n")

    for lnum in limited_num:
        
        f.write("curl -O ftp://stripero:s0l0le0@ftp.iac.es/coadds/f"+lnum+"_u.rec.fits\n")
        f.write("gzip f"+lnum+"_u.rec.fits\n")
        f.write("curl -O ftp://stripero:s0l0le0@ftp.iac.es/coadds/f"+lnum+"_g.rec.fits\n")
        f.write("gzip f"+lnum+"_g.rec.fits\n")
        f.write("curl -O ftp://stripero:s0l0le0@ftp.iac.es/coadds/f"+lnum+"_r.rec.fits\n")
        f.write("gzip f"+lnum+"_r.rec.fits\n")
        f.write("curl -O ftp://stripero:s0l0le0@ftp.iac.es/coadds/f"+lnum+"_i.rec.fits\n")
        f.write("gzip f"+lnum+"_i.rec.fits\n")
        f.write("curl -O ftp://stripero:s0l0le0@ftp.iac.es/coadds/f"+lnum+"_z.rec.fits\n")
        f.write("gzip f"+lnum+"_z.rec.fits\n")

f.close()

