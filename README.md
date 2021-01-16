# MPAS-SW-TL-AD

> NOTE: This library is still under active development, any suggestions and comments are greatly appreciated. It has been extensively tested with idealized and real atmosphere (ERA5 500 hPa) simulations.


## Prerequisites
---
* [Anaconda](https://www.anaconda.com/) (recommended)
* [F90Nml](https://pypi.org/project/f90nml/)
* [netCDF4](https://pypi.org/project/netCDF4/)

This package includes the MPAS-Shallow Water dynamics under the Python-Fortran structure and its tangent linear and adjoint components [1].

## Installation
---
Run the following commands to compile the package:
    cd src
    make
Successful compilation should have a module_sw_mpas.so or module_sw_mpas.*.so generated for future Python import.

[1]: Tian, X. (2020). Evolutions of Errors in the Global Multiresolution Model for Prediction Across Scales - Shallow Water (MPAS-SW). *Q. J. Royal Meteorol. Soc.*, [https://doi.org/10.1002/qj.3923](https://doi.org/10.1002/qj.3923)
