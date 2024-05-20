# Stripe82 cutouts & mocks
Building galaxy segmentation imaging datasets from Stripe82 SDSS. This repository includes the code to produce cutouts from the observations plus their segmentation masks and to make stamps of mock galaxies simulated as bulge+disc light profiles.

The segmentation masks are delimited by R_edge, radius at which the galaxy has a truncation in its disc (given by a change in the slope of its radial light profile).

## Stripe82 cutouts
The galaxies considered are those studied in [Chamba+ 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...667A..87C/abstract). You'll need access to the catalogues produced for this article, which can be found [here](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/667/A87) (Dwarfs catalogue and Early-/Late-types catalogue).

The observation frames can be downloaded from the IAC [Stripe82 website](http://research.iac.es/proyecto/stripe82/). You can use the script *** in this repository to dowload only those frames that include the galaxies from a given catalogue. Additional frames may be needed depending of the cutout size chosen. This can be solved by iterating a few times the codes *** and ***.

### Steps to produce the cutouts


## Stripe82 mock dataset
The mock galaxies consist of a combination of bulge + disc profiles produced with [Imfit](https://www.mpe.mpg.de/~erwin/code/imfit/) (required) and convolved with the Sloan PSFs in the u,g,r,i and z bands. These are later inserted in a real Stripe82 image.

### Steps to produce the mocks
[PSFs](https://ui.adsabs.harvard.edu/abs/2020RNAAS...4..130I/abstract)


