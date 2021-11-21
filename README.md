# Helical and skyrmion dynamics in MnSi


## Description
This is the supplementary source code to the paper
*Topological magnon band structure of emergent Landau levels in a skyrmion lattice*.

The code contains a *Takin* plugin module and its helper tools.
The module calculates the dispersion and the dynamical structure factor for the
helimagnetic, the field-polarised, and the skyrmion phase of MnSi.

The development repository can be found here:
https://code.ill.fr/scientific-software/takin/plugins/mnsi


## Dependencies
[*Takin 2*](https://doi.org/10.5281/zenodo.4117437) and the *tlibs*/*tlibs2* libraries are needed for compilation.
Their repositories are available here:
	- https://github.com/t-weber/takin2
	- https://code.ill.fr/scientific-software/takin


## Setup
	- Download external dependencies: `cd ext && ./setup_externals.sh && cd ..`.
	- The ext/ directory contains the source code of the external libraries.
	- Build the module: `make -j4`.


## Acknowledgements
This code is based on theoretical magnon dispersion models and their Mathematica
implementations by M. Garst and J. Waizner, see our papers above and these references:
	- https://doi.org/10.1088/1361-6463/aa7573
	- https://kups.ub.uni-koeln.de/7937/

The helimagnon and ferromagnetic part of this code has been used in the following papers:
	- https://doi.org/10.1103/PhysRevB.100.060404
	- https://doi.org/10.1103/PhysRevB.97.224403
	- https://doi.org/10.1063/1.5041036

Furthermore, this code is based on Python implementations by M. Kugler and G. Brandl of
early versions of the mentioned theoretical models.
