#
# calculates high symmetry points for the skx lattice
# @author Tobias Weber <tweber@ill.fr>
# @date 9-feb-2023
# @license GPLv2 (see 'LICENSE' file)
#

import LinearAlgebra as la
using Printf


# field, pinning and lattice vector
Bvec = [ 1;  1; 0 ]
Pvec = [ 1; -1; 0 ]
Gvec = [ 1;  1; 0 ]
include_G = false


#
# Rodrigues' formula to rotate vec1 into vec2
# @see https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
# @see (Merziger 1993), p. 208
# @see (Arens 2015), p. 718 and p. 816
#
function rotmat(vec1, vec2)
	axis = la.cross(vec1, vec2)
	axislen = la.norm2(axis)
	axis /= axislen
	angle = atan(axislen, la.dot(vec1, vec2))
	#@printf("axis = %s, angle = %f\n", axis, angle/pi*180)

	s = sin(angle)
	c = cos(angle)

	skew = [
		       0  -axis[3]   axis[2];
		 axis[3]         0  -axis[1];
		-axis[2]   axis[1]         0 ]

	M = (1 - c) * axis*la.transpose(axis)
	M += la.diagm([c; c; c]) + s*skew
	return M
end


Brot = rotmat(Bvec, [0; 0; 1])
Pvec = Brot * Pvec
Prot = rotmat(Pvec, [1; 0; 0])

rot = Prot * Brot
rotinv = la.inv(rot)
#@printf("rot = %s\ninv_rot = %s\n", rot, rotinv)


B = [1 cos(120/180*pi) 0; 0 sin(120/180*pi) 0; 0 0 1]

posK1 = B * [1; 0; 0]
posK1 /= la.norm2(posK1)

posK2 = B * [1; 1; 0]
posK2 /= la.norm2(posK2)

posM = (posK1 + posK2) / 2

posK1 = rotinv * posK1
posK2 = rotinv * posK2
posM = rotinv * posM

posK1 *= 0.039 / (2*pi / 4.558)
posK2 *= 0.039 / (2*pi / 4.558)
posM *= 0.039 / (2*pi / 4.558)


if include_G
	posK1 += Gvec
	posK2 += Gvec
	posM += Gvec
end


@printf("K1 = %s\n", round.(posK1, digits=5))
@printf("K2 = %s\n", round.(posK2, digits=5))
@printf("M  = %s\n", round.(posM, digits=5))
