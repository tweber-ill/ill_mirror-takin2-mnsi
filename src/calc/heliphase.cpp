/**
 * Calculates the helimagnetic phase diagram
 * @author tweber@ill.fr
 * @date aug-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/heli.h"


int main()
{
	using t_real = double;
	using t_cplx = std::complex<t_real>;
	const auto j = t_cplx(0,1);

	Heli<t_real, t_cplx, 4> heli;
	std::vector<ublas::vector<t_cplx>> fourier{
		tl2::make_vec<ublas::vector<t_cplx>>({0, 0, 0.1}),
		// helical order => Re{M} perp. Im{M}
		tl2::make_vec<ublas::vector<t_cplx>>({1.+j, 1.-j, 0}) / std::sqrt(2),
	};

	heli.SetFourier(fourier);
	heli.SetT(-100);
	heli.SetB(0);
	heli.SaveStates("heli.dat", 4, 0, 0, 0, 0);

	return 0;
}
