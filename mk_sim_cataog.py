"""Things to change:
    - when the height of the galaxy (z0) is big relative to the Redge (over 10%)
    imfit produces a disc-like thing going in the perpendicular direction.
    
    - limit the redge considering the redshift, such that the distribution of z_gal
    and redge (in kpc) remains uniform."""
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.cosmology.realizations import Planck15 as cosmo
import numpy as np
import matplotlib.pyplot as plt
import yaml

with open("config.yaml", "r") as f: 
    config = yaml.load_all(f)

rng = np.random.default_rng(12345)
# main properties of the galaxy


n_gal = 1000
pix_scale = 0.396 * u.arcsec
z_min = 0.003
z_max = 0.09
img_size = config["sim"]["cutout_size"]
pix_margin = 20


rmax = img_size/2.-pix_margin
zs = np.linspace(z_min, z_max, 100000)
rmax_kpc = rmax * pix_scale / cosmo.arcsec_per_kpc_proper(zs)
redshift_bottom = zs[np.searchsorted(rmax_kpc, 70.*u.kpc)-1]
redshift_upper = z_max-redshift_bottom
plt.plot([z_min, z_max], [1.5,1.5])
plt.plot([z_min, z_max], [70,70])
plt.plot(zs,rmax_kpc)
plt.plot(zs+z_max-redshift_bottom,rmax_kpc.value- (rmax * pix_scale.value / cosmo.arcsec_per_kpc_proper(z_min).value))
plt.show()



z_gal = rng.uniform(z_min, z_max, size=n_gal)
fraction_b = rng.uniform(0.01, 0.2, size=n_gal)

# truncated disc parameters
pa_disc = rng.uniform(0., 180., size=n_gal) #* u.deg
ell_disc = rng.uniform(0., 0.7, size=n_gal) # 1- (b/a) or 1-q_ellip
incl_disc = rng.uniform(0., 90., size=n_gal)# * u.deg #180./np.pi * np.arccos(1. - ell_)
r_edge = rng.uniform(1.5, 70., size=n_gal) #* u.kpc # in kpc
r_edge = np.where(z_gal<redshift_bottom, 
                  rng.uniform(1.5, rmax * pix_scale.value / cosmo.arcsec_per_kpc_proper(z_gal).value, size=n_gal),
                  r_edge)
r_min = rmax * pix_scale.value / cosmo.arcsec_per_kpc_proper(z_gal-redshift_upper).value \
        - (rmax * pix_scale.value / cosmo.arcsec_per_kpc_proper(z_min).value)
r_min = np.where(r_min<0, 0, r_min)
r_edge = np.where(z_gal>redshift_upper, 
                  rng.uniform(r_min, 70., size=n_gal),
                  r_edge)
h1_ = r_edge * rng.uniform(0.2, 0.4, size=n_gal)
h2_ = h1_ * rng.uniform(0.25, 0.35, size=n_gal)
alpha = rng.uniform(0.5, 7.5, size=n_gal) #/ u.kpc# in 1/kpc
z0 = np.array([rng.uniform(0.05*r_edge[ii],0.08*r_edge[ii]) for ii in range(n_gal)])#(0.7,1.5, size=n_gal) #* u.kpc# in kpc
n_disc =  np.ones(n_gal, dtype=np.single)

plt.scatter(z_gal, r_edge)
plt.show()
# bulge parameters
PA_bulge = 0 * np.ones(n_gal, dtype=np.single)
ell_bulge = 0 * np.ones(n_gal, dtype=np.single)
Re_bulge = rng.uniform(0.5, 4., size=n_gal)#np.array([rng.uniform(0.5, 0.2*rr) for rr in r_edge]) * u.kpc
n_bulge = rng.uniform(1., 5., size=n_gal)#4 * np.ones(n_gal, dtype=np.single)


# Assign observed sed
hdu_obs = fits.open(config["obs"]["catalogue"])
gal_id = hdu_obs[1].data["JID"]
sed_ind = rng.choice(gal_id.shape[0], size=n_gal)
sed_factor = rng.uniform(0.9, 1.1, size=n_gal)
obs_sed = [gal_id[ii] for ii in sed_ind]
f_seds = np.load(config["sim"]["seds_path"]+"seds_obs.npy")
fseds_assign = np.tile(sed_factor, (5,1)).T * f_seds[sed_ind]


gal_prop_t = Table([pa_disc, ell_disc, incl_disc, r_edge, h1_, h2_,
                     alpha, z0, n_disc, PA_bulge, ell_bulge, Re_bulge,
                     n_bulge, z_gal, obs_sed, sed_ind, sed_factor, fraction_b, 
                     fseds_assign], 
                   names=["pa_disc", "ell_disc", "incl_disc", "r_edge", "h1_", 
                    "h2_", "alpha", "z0", "n_disc", "PA_bulge", "ell_bulge",
                    "Re_bulge", "n_bulge", "z_gal", "obs_sed", "sed_ind", 
                    "sed_factor", "fraction_b", "target_sed"])



hdu = fits.BinTableHDU(data=gal_prop_t)
hdu.writeto(config["sim"]["catalogue"], overwrite=True)  
