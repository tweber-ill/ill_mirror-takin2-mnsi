#
# demagnetisation factors
# @author Tobias Weber <tweber@ill.fr>
# @date apr-2022
# @license GPLv2 (see 'LICENSE' file)
#

using Printf

# -----------------------------------------------------------------------------
# constants
# The values are from the theoretical models by. M. Garst and J. Waizner:
#      - https://doi.org/10.1088/1361-6463/aa7573
#      - https://kups.ub.uni-koeln.de/7937/
#      - Personal communications with M. Garst, 2017-2020.
# Further values, e.g. hoc and chi, are from the paper and its Python code:
#      - https://doi.org/10.1103/PhysRevLett.115.097203
# -----------------------------------------------------------------------------
chi = 0.34
# -----------------------------------------------------------------------------


#
# cylindrical demagnetisation factor
# @see https://doi.org/10.1063/1.343481
#
function demag_cyl_x(l, r)
	sqrt_area = sqrt(pi)*r
	Nx = l / (2*l + sqrt_area)
	return Nx
end


#
# cylindrical demagnetisation factor
# @see https://doi.org/10.1063/1.343481
#
function demag_cyl_z(l, r)
	sqrt_area = sqrt(pi)*r
	Nz = sqrt_area/(2*l + sqrt_area)
	return Nz
end


l = 30
r = 5

Nx = demag_cyl_x(l, r)
Nz = demag_cyl_z(l, r)

chi_ratio_x = 1 / (1 + Nx * chi)
chi_ratio_z = 1 / (1 + Nz * chi)

@printf("Nx = %.4f, Nz = %.4f\n", Nx, Nz)
@printf("rx = %.4f, rz = %.4f\n", chi_ratio_x, chi_ratio_z)
