# RadioLensfit2

**Radio Weak Lensing shear measurement in the Visibility Domain**

This is version 2.0 of [RadioLensfit](https://github.com/marziarivi/RadioLensfit).
Three tools are provided:

**Simulate**

Simulation of a radio weak lensing observation with a given reduced cosmic shear and a source flux threshold Fmin.
Radio telescope configuration and observing time sampling must be provided in a Measurement Set.

The simulated faint radio catalog is populated by star-forming galaxies according to flux, scalelength, ellipticity distributions estimated by JVLA observations (SWIRE and COSMOS). The number of galaxies is dependent on the flux prior and provided flux-cut.
Sources are simulated according the RING TEST to avoid shape noise: same source flux and size for a fixed number of ellipticity values simmetrically distributed on the same ring. The effect of the input reduced cosmic shear is applied to galaxies ellipticity.
This source catalog is written in a text file called *RWL_galaxy_catalog_<FoV>_<Fmin>.txt*, where FoV is the provided effective field of view. 
The corresponding simulated visibilities observed by the radio telescope are written in the DATA column of the same Measurement Set (I stokes component). The instrument noise is added. NO primary beam effect is currently considered.
  
Usage: *Simulate* (filename MS) (effective field of view [arcmin]) (min flux [muJy]) (shear coord_1) (shear coord_2)
 
**RadioLensfit2-MS**

Measurement of star forming galaxy ellipticies from a radio weak lensing observation.
Data visibilities and observation configuration must be provided in a Measurement Set. 
The number of galaxies and the corresponding source catalog (ordered by decreasing flux) containing source SNR, position and flux must be provided. Source position and flux is used for source visibilities extraction, then shape measurement is performed according to RadioLensfit methodology: a single model fitting approach where the likelihood is marginalised over position, flux and scalelength source parameters.  
The list of galaxies with the measured ellipticities is written in a text file called *ellipticities.txt*.

Usage: *RadioLensfit2-MS* (filename MS) (filename source catalog) (num_sources)
 
**shear.py** 

Shear measurement from weighted average of the measured ellipticity of the galaxies in the field of view.

Usage: *shear.py* (filename measured ellipticities)

# Installation

Requirements:
- CASACORE library for Measurement Set I/O
- GSL library and 

1. Edit the Makefile:

enable/disable OpenMP (default: enabled)

enable/disable gridding (default: gridding enabled)

set the compiler and compilation flags you want to use (default: GNU)

2. make all

# Citing RadioLensfit2

If you use RadioLensfit2 and find it useful, please consider citing the related papers:

Rivi M., Miller L., Makhathini S., Abdalla F. B., 2016, MNRAS, 463, 1881 - [arXiv:1603.04784](https://arxiv.org/abs/1603.04784)

Rivi M., Miller L., 2018, MNRAS, 476, 2053 - [arXiv:1709.01827](https://arxiv.org/abs/1709.01827)
