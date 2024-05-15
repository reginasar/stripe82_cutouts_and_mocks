# Stripe82 cutouts & mocks
Building galaxy segmentation imaging datasets from Stripe82 SDSS. This repository includes the code to produce cutouts from the observations plus their segmentation masks and to make stamps of mock galaxies simulated as bulge+disc light profiles.

The segmentation masks are delimited by R_edge, radius at which the galaxy has a truncation in its disc (given by a change in the slope of its radial light profile).

## Stripe82 cutouts
The galaxies considered are those studied in Chamba+ 2022. You'll need access to the catalogues produced for this article, which can be found here.

The observation frames can be downloaded from the IAC website. You can use the script *** in this repository to dowload only those frames that include the galaxies from a given catalogue. Additional frames may be needed depending of the cutout size chosen. This can be solved by iterating a few times the codes *** and ***.

### Steps to produce the cutouts


## Stripe82 mock dataset
The mock galaxies consist of a combination of bulge + disc profiles produced with imfit (required) and convolved with the Sloan PSFs in the u,g,r,i and z bands. These are later inserted in a real Stripe82 image.

### Steps to produce the mocks



