# RadioLensfit2

**Radio Weak Lensing shear measurement in the Visibility Domain**

This is version 2.0 of [RadioLensfit](https://github.com/marziarivi/RadioLensfit).
Three tools are provided:

**Simulate**

Simulation of a radio weak lensing observation with a given reduced cosmic shear and a source flux threshold Fmin (or a given source catalog).
Radio telescope configuration and observing time sampling must be provided in a Measurement Set.

The simulated faint radio catalog is populated by star-forming galaxies according to flux, scalelength, ellipticity distributions estimated by JVLA observations (SWIRE and COSMOS). The number of galaxies is dependent on the flux prior and provided flux-cut.
Sources are simulated according to the RING TEST to avoid shape noise: same source flux and size for a fixed number of ellipticity values simmetrically distributed on the same ring. The effect of the input reduced cosmic shear is applied to galaxies ellipticity.
This source catalog is written in a text file called *RWL_galaxy_catalog_(FoV)_(Fmin).txt*, where FoV is the provided effective field of view. 
The corresponding simulated visibilities observed by the radio telescope are written in the DATA column of the same Measurement Set (I stokes component). The instrument noise is added. NO primary beam effect is currently considered.
  
Usage: *Simulate* (filename MS) (effective field of view [arcmin]) (min flux [muJy]) (shear coord_1) (shear coord_2)
Usage: *Simulate* (filename MS) (source catalog filename) (number of sources) (shear coord_1) (shear coord_2)

The version of simulate must be selected in the Makefile
 
**RadioLensfit2**

Measurement of star forming galaxy ellipticies from a radio weak lensing observation.
Data visibilities and observation configuration must be provided in a Measurement Set. 
The number of galaxies and the corresponding source catalog (ordered by decreasing flux) containing source SNR, position and flux must be provided. Source position and flux is used for source visibilities extraction, then shape measurement is performed according to RadioLensfit methodology: a single model fitting approach where the likelihood is marginalised over position, flux and scalelength source parameters.  
The list of galaxies with the measured ellipticities is written in a text file called *ellipticities.txt*.

Serial version usage: *RadioLensfit2* (source catalog filename) (number of sources) (filename MS).MS

MPI version usage: *RadioLensfit2-mpi* (source catalog filename) (number of sources) (filename MSs prefix).MS

MS must be split in individual spectral windows. All MS must have the same name ending with the number of the spectral window.
Filename prefix consists in the common part of all the MS name, i.e. (without the final number and ".MS"   
 
**shear.py** 

Shear measurement from weighted average of the measured ellipticity of the galaxies in the field of view. 
To be modified...


# Installation

Requirements:
- CASACORE library for Measurement Set I/O
- GSL library 
- MPI library (for the parallel version)

Possibly you can modify the fixed parameters defined in the *default_params.h* file

Edit the Makefile:

- enable/disable OpenMP (default: enabled)
- enable/disable MPI (default: enabled)

- update CASACORE include and library path
- set the compiler and compilation flags you want to use (default: GNU)

*make all*

# Citing RadioLensfit2

If you use RadioLensfit2 and find it useful, please consider citing the related papers:

Rivi M., Miller L., Makhathini S., Abdalla F. B., 2016, MNRAS, 463, 1881 - [arXiv:1603.04784](https://arxiv.org/abs/1603.04784)

Rivi M., Miller L., 2018, MNRAS, 476, 2053 - [arXiv:1709.01827](https://arxiv.org/abs/1709.01827)
