![alt text](https://github.com/dotbot2000/iso/blob/master/plots/iso.png)

## Overview

This project provides a means to transform [MESA](http://mesa.sourceforge.net) history files into a uniform basis for interpolation and then construct new stellar evolution tracks and isochrones from that basis. It is developed in Fortran.

The theory is described in [this ApJS paper](http://adsabs.harvard.edu/abs/2016ApJS..222....8D).  You are encouraged to read the paper before attempting to use these programs. If you use these programs in your research, please cite the paper!

## Prerequisites

(Note: If you already have MESA installed, then you most likely meet all these requirements!)

1. `git` if you want to clone the repository.  Alternatively, you can download a zip file from the github page.

2. `make` to handle compilation.  make is distributed as a package in Linux distributions and comes as part of the Xcode Command Line Tools in Mac OS X.

3. The code uses MESA libraries for things like interpolation so you'll need to have MESA installed--at least the utilities and numerical libraries.  See [here](http://mesa.sourceforge.net/prereqs.html) for instructions.

4. A modern Fortran compiler (but if you have MESA installed already, you'll have it).

5. In order to identify primary EEPs in MESA history files, several history columns are required:
   + `star_age`
   + `star_mass`
   + `log_LH`
   + `log_LHe` 
   + `log_Teff` 
   + `log_L`
   + `log_g` 
   + `log_center_T` 
   + `log_center_Rho` 
   + `center_h1`
   + `center_he4`
   + `center_c12`
   + `center_gamma`
   + `surface_h1`
   + `he_core_mass`
   + `c_core_mass`

6. The `history_columns.list` file that was used in generating the MESA history files. By default it resides in `$MESA_DIR/star/defaults/history_columns.list`.

7. Synthetic color-magnitude diagrams require a set of bolometric correction tables. These will be made available coincident with the release of MIST Paper I.

## Installation

It's a simple matter of cloning the repository (for which you'll need git installed) or downloading a zip file from github.  

Then cd into the project directory and run the mk script to compile the codes.

```
git clone git@github.com:dotbot2000/iso.git

cd iso

./mk
```

## Getting Started

The basic workflow is split into two parts:

1. Convert the MESA history files into EEP-based tracks (.eep files) via the program `make_eep`.

2. Use the EEP files to create new stellar evolution tracks with `make_track` and/or isochrones with `make_iso`.

Both programs read the same input file. An input file is divided into two sections, a track/eep section and an ioschrone section.  The track/eep section comes first; it contains:
+ a list of data file directories to which the program will both read from and write to.  
+ the MESA history columns file "history_columns.list" that was used to create the MESA history files (tracks)
+ the number of tracks
+ a list of the tracks sorted by increasing order of initial mass.

The isochrone section follows the track/eep section; it contains:
+ the isochrone output file name
+ one of two options that specify how the ages will be input, either `min_max` or `list`.
  If `min_max` then you'll give the desired number of isochrones, the minimum age and the maximum age on three successive lines.
  If `list` then you'll give the desired number of isochrones followed by a list of individual ages, one per line.
+ The scale of ages, either `linear` or `log10`, in which the isochrone ages will be specified. Isochrone ages are always given in years. For example, 10 Gyr could be specified as either 1e10 (linear) or 10.0 (log10).

See the file `input.example` for the general layout.

A variety of code options can be set via two namelists, both reside in the file `input.nml`. The number of secondary EEPs that are desired between each pair of primary EEPs is set in the file `input.eep`. If this file does not exist, the code will set a default (small) number of 50 secondary EEPs between each pair of primary EEPs.


After these input files are configured correctly, run the codes using the following commands:

`./make_eep input.example`

check terminal output for error messages.

`./make_iso input.example`

If all goes well, then you'll have a fresh set of isochrones to explore.

### NOTES 

* Sometime around MESA revision 5000, the name of the "core mass" history columns changed from, e.g., `h1_boundary_mass` (old) to `he_core_mass` (new). This is dealt with in the code by the variable `old_core_mass_names`: to use the old names, set this to `.true.`; otherwise leave it `.false.`.
