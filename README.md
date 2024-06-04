# Stripe82 cutouts & mocks
Building galaxy segmentation imaging datasets from Stripe82 SDSS. This repository includes the code to produce cutouts from the observations plus their segmentation masks and to make stamps of mock galaxies simulated as bulge+disc light profiles.

The segmentation masks are delimited by R_edge, radius at which the galaxy has a truncation in its disc (given by a change in the slope of its radial light profile).

To make the steps easier, set up your own config file (see instructions at the end of repository). 

## Stripe82 cutouts
The galaxies considered are those studied in [Chamba+ 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...667A..87C/abstract). You'll need access to the catalogues produced for this article, which can be found [here](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/667/A87) (Dwarfs catalogue and Early-/Late-types catalogue).

The observation frames can be downloaded from the IAC [Stripe82 website](http://research.iac.es/proyecto/stripe82/). Go to step 1 in this section to dowload only those frames that include the galaxies from a given catalogue. Additional frames may be needed depending of the cutout size chosen (step 2).


### Steps to produce the cutouts

1) Download the frames required to produce the cutouts. Run:

python mk_download_file.py
bash download.sh

2) Produce the cutout stamps witht their respective masks. Additional frames may be needed depending of the cutout size chosen. This can be solved by running download_update.sh and re-running the cutout generator. Run:

python mk_stripe82_stamps_beta.py
bash download_update.sh
python mk_stripe82_stamps_beta.py

At this point your observed galaxy stamps are ready.

## Stripe82 mock dataset
The mock galaxies consist of a combination of bulge + disc profiles produced with [Imfit](https://www.mpe.mpg.de/~erwin/code/imfit/) (required) and convolved with the Sloan PSFs in the u,g,r,i and z bands. These are later inserted in a real Stripe82 image.

### Steps to produce the mocks

1) Subtract realistic SEDs from the observed galaxies with:

'''
python get_seds.py
'''

2) Make a catalogue with properties of your simulated galaxies. Run:

'''
python mk_catalog.py
'''

3) Get your PSFs. Either with a simple gaussian approximation or downloading realistic SDSS PSFs from here:

[PSFs](https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.5317I/abstract)

Make sure that the files do not contain NaNs (replace with zeros if necessary). To accelarete the convolution step, reshape your PSFs to have dimensions smaller or equal to '2N+1', where N is the side of the stamps you want to make.

4) Make the simulated profiles for all the galaxies defined in the simulated catalog (step 2). Run:

'''
python mk_imfit_sim_beta.py
'''

5) Get back/foreground frames. In case you don't want to re-use the frames of the observed galaxies, you should download a new set of frames with: 


Otherwise skip this step and just set the bg_frames option in the config file to be the path of the observed frames used earlier. 

6) Add realism to your simulated stamps (step 4) by adjusting the flux of the simulated profiles to realistic values (obtained with step 1) and overlap with an observed back/foreground (previous step). Run:

'''
python mk_mag_fits_beta.py
'''


## Setup config file

General guidelines:

- If the key indicates a path make sure your value is a string with an existent path. The string must finish with the "/" character.

- There are two broad sections "obs" and "sim". Options in "sim" will not be used to produce stamps of observed galaxies. However, some "obs" options wil be used to generate mock stamps. 

- The config file should be completed as:

'''
obs:
  catalogue: "/path/to/directory/catalogue_obs_galaxies.fits"
  download_script: "/path/to/directory/download.sh"
  stripe82_frames_path: "/path/to/directory/with/frames/"
  cutout_size: 600
  out_path: "/path/to/directory/with/cutouts/"
  seds_path: "/path/to/directory/"

sim:
  catalogue: "/path/to/directory/catalogue_sim_galaxies.fits"
  config_files_path: "/path/to/directory/with/config_files/for/simulations/"
  simulated_files_path: "/path/to/directory/with/simulated/profiles/stamps/"
  psf_file: "/path/to/directory/psf-filter-?-v0.5.fits[1]"
  make_image: "/path/to/diectory/makeimage"
  background_frames_path: "/path/to/directory/bg_frames/"
  n_gal: 1000
  cutout_size: 600
  gen_file: "/path/to/directory/gen_fits.sh"
  out_path: "/path/to/directory/mocks/"

'''
  

- "?" in psf_file will be replaced with the band name (u,g,r,i or z). The number in "[]" indicates the extensionwhere the PSF is, if Primary then remove "[]".



