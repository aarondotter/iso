![alt text](https://github.com/dotbot2000/iso/blob/master/plots/iso.png)

## Overview

This project provides a means to transform [MESA](http://mesa.sourceforge.net) history files into a uniform basis for interpolation and then construct new stellar evolution tracks and isochrones from that basis. It is developed in Fortran.

The theory is described in [this paper]().

## Prerequisites

1. `git` if you want to clone the repository.  Alternatively, you can download a zip file from the github page.

2. `make` to handle compilation.  make is distributed as a package in Linux distributions and comes as part of the Xcode Command Line Tools in Mac OS X.

3. The code uses MESA libraries for things like interpolation so you'll need to have MESA installed--at least the utilities and numerical libraries.  See [here](http://mesa.sourceforge.net/prereqs.html) for instructions.

4. A modern Fortran compiler (but if you have MESA installed already, you'll have it).

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


