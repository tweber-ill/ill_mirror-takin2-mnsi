# Skyrmion, helical, and field-polarised magnon dynamics in MnSi

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5718363.svg)](https://doi.org/10.5281/zenodo.5718363)


## Description
This is the supplementary source code to the paper
*Topological magnon band structure of emergent Landau levels in a skyrmion lattice*.

The code contains a [*Takin*](https://doi.org/10.5281/zenodo.4117437) plugin module and its helper tools.
The module calculates the dispersion and the dynamical structure factor for the helimagnetic, the field-polarised, and the skyrmion phase of MnSi.

The development repository can be found here:
- https://code.ill.fr/scientific-software/takin/plugins/mnsi


## Dependencies
The [*Takin* software](https://doi.org/10.5281/zenodo.4117437) and the [*tlibs* libraries](https://doi.org/10.5281/zenodo.5717779) are needed for compilation.
Their repositories are available here:
- Stable releases: https://github.com/t-weber/takin2
- Development versions: https://code.ill.fr/scientific-software/takin
- Binary releases: https://wiki.mlz-garching.de/takin
- DOIs: Takin 2: [10.5281/zenodo.4117437](https://doi.org/10.5281/zenodo.4117437), tlibs 2: [10.5281/zenodo.5717779](https://doi.org/10.5281/zenodo.5717779), old Takin 1 and tlibs 1: [10.5281/zenodo.3961491](https://doi.org/10.5281/zenodo.3961491).


## Setup
- Download the external dependencies: `cd ext && ./setup_externals.sh && cd ..`.
- The ext/ directory should now contain the source code of the external libraries.
- Build the module: `make -j4`.
- Copy the built plugin modules to *Takin's* plugin directory: `mkdir -pv ~/.takin/plugins/ && cp -v lib/*.so ~/.takin/plugins/`
- The helper tools can be found in the bin/ directory.


## Acknowledgements
This source code is based on theoretical magnon dispersion models and their *Mathematica* implementations by M. Garst and J. Waizner, see these references and our papers below:
- M. Garst and J. Waizner, Skyrmion linear spin-wave theory and *Mathematica* implementation, personal communications (2017-2020).
- M. Garst and J. Waizner, Helimagnon linear spin-wave model and *Mathematica* implementation, personal communications (2014-2019).
- M. Garst and J. Waizner, Field-polarised linear spin-wave model and *Mathematica*, personal communications (2016-2019).
- M. Garst, J. Waizner, and D. Grundler, J. Phys. D: Appl. Phys. **50** 293002, https://doi.org/10.1088/1361-6463/aa7573 (2017).
- J. Waizner, PhD thesis, Universität zu Köln, https://kups.ub.uni-koeln.de/7937/ (2017).

Furthermore, the present source code is based on optimised *Python* implementations by M. Kugler and G. Brandl of early versions of the theoretical models mentioned above; it started as a translation of the following *Python* codes into *C++*:
- G. Brandl and M. Kugler, Helimagnon implementation in *Python*, personal communications (2015-2016).
- M. Kugler, G. Brandl, J. Waizner, M. Janoschek, R. Georgii, A. Bauer, K. Seemann, A. Rosch, C. Pfleiderer, P. Böni, and M. Garst, Phys. Rev. Lett. **115**, 097203, https://doi.org/10.1103/PhysRevLett.115.097203 (2015).
- M. Kugler and G. Brandl, Skyrmion spin-wave implementation in *Python*, personal communication (2016).

The following alternate *Python* implementations of the skyrmion spin-wave model exist:
- D. Fobes, Implementation in *Python*, personal communication (2016).
- L. Beddrich, Implementation in *Python*, https://github.com/LukasBeddrich/skyrmion-model (2017).

The helimagnon and ferromagnetic parts of this code have been used in the following papers:
- T. Weber, J. Waizner, P. Steffens, A. Bauer, C. Pfleiderer, M. Garst, and P. Böni, Phys. Rev. B **100**, 060404(R), https://doi.org/10.1103/PhysRevB.100.060404 (2019).
- T. Weber, J. Waizner, G. S. Tucker, R. Georgii, M. Kugler, A. Bauer, C. Pfleiderer, M. Garst, and P. Böni, Phys. Rev. B **97**, 224403, https://doi.org/10.1103/PhysRevB.97.224403 (2018).
- T. Weber, J. Waizner, G. S. Tucker, L. Beddrich, M. Skoulatos, R. Georgii, A. Bauer, C. Pfleiderer, M. Garst, and P. Böni, AIP Advances **8**, 101328, https://doi.org/10.1063/1.5041036 (2018).
