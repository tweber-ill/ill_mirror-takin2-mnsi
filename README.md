# Skyrmion and helical dynamics in MnSi


## Description
This is the supplementary source code to the paper
*Topological magnon band structure of emergent Landau levels in a skyrmion lattice*.

The code contains a *Takin* plugin module and its helper tools.
The module calculates the dispersion and the dynamical structure factor for the helimagnetic, the field-polarised, and the skyrmion phase of MnSi.

The development repository can be found here:
- https://code.ill.fr/scientific-software/takin/plugins/mnsi


## Dependencies
[*Takin 2*](https://doi.org/10.5281/zenodo.4117437) and the *tlibs* / *tlibs2* libraries are needed for compilation.
Their repositories are available here:
- https://github.com/t-weber/takin2
- https://code.ill.fr/scientific-software/takin


## Setup
- Download external dependencies: `cd ext && ./setup_externals.sh && cd ..`.
- The ext/ directory contains the source code of the external libraries.
- Build the module: `make -j4`.


## Acknowledgements
This code is based on theoretical magnon dispersion models and their *Mathematica* implementations by M. Garst and J. Waizner, see these references and our papers below:
- M. Garst and J. Waizner, Skyrmion linear spin-wave theory and *Mathematica* implementation, personal communications, 2017-2020.
- M. Garst and J. Waizner, Helimagnon linear spin-wave model and *Mathematica*, personal communications, 2014-2018.
- M. Garst and J. Waizner, Field-polarised linear spin-wave model and *Mathematica*, personal communications, 2016.
- M. Garst, J. Waizner, and D. Grundler, J. Phys. D: Appl. Phys. **50** 293002, https://doi.org/10.1088/1361-6463/aa7573 (2017).
- J. Waizner, PhD thesis, Universität zu Köln, https://kups.ub.uni-koeln.de/7937/ (2017).

Furthermore, the present code is based on *Python* implementations by M. Kugler and G. Brandl of early versions of the theoretical models mentioned above:
- G. Brandl and M. Kugler, Helimagnon implementation in *Python*, personal communication (2015-2016)
- M. Kugler, G. Brandl, J. Waizner, M. Janoschek, R. Georgii, A. Bauer, K. Seemann, A. Rosch, C. Pfleiderer, P. Böni, and M. Garst, Phys. Rev. Lett. **115**, 097203, https://doi.org/10.1103/PhysRevLett.115.097203 (2015).
- M. Kugler and G. Brandl, Skyrmion spin-wave implementation in *Python*, personal communication (2016).

The following alternate *Python* implementations of the skyrmion spin-wave model exist:
- D. Fobes, Implementation in *Python*, personal communication (2016).
- L. Beddrich, Implementation in *Python*, https://github.com/LukasBeddrich/skyrmion-model (2017).

The helimagnon and ferromagnetic parts of this code have been used in the following papers:
- T. Weber, J. Waizner, P. Steffens, A. Bauer, C. Pfleiderer, M. Garst, and P. Böni, Phys. Rev. B **100**, 060404(R), https://doi.org/10.1103/PhysRevB.100.060404 (2019).
- T. Weber, J. Waizner, G. S. Tucker, R. Georgii, M. Kugler, A. Bauer, C. Pfleiderer, M. Garst, and P. Böni, Phys. Rev. B **97**, 224403, https://doi.org/10.1103/PhysRevB.97.224403 (2018).
- T. Weber, J. Waizner, G. S. Tucker, L. Beddrich, M. Skoulatos, R. Georgii, A. Bauer, C. Pfleiderer, M. Garst, and P. Böni, AIP Advances **8**, 101328, https://doi.org/10.1063/1.5041036 (2018).
