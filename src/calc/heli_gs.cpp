/**
 * Calculate the ground state fourier components and free energy
 * @author tweber@ill.fr
 * @date aug-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/heli.h"
#include <fstream>

#define ORDER DEF_HELI_ORDER

using t_real = double;
using t_cplx = std::complex<t_real>;


int main()
{
	const auto j = t_cplx(0, 1);

	Heli<t_real, t_cplx, DEF_HELI_ORDER> heli;
	std::vector<ublas::vector<t_cplx>> fourier{
		tl2::make_vec<ublas::vector<t_cplx>>({0, 0, 1.}),
		// helical order => Re{M} perp. Im{M}
		tl2::make_vec<ublas::vector<t_cplx>>({1.+j, 1.-j, 0}) / std::sqrt(2),
 	};

        heli.SetFourier(fourier);
        heli.SetT(-1000, false);
        heli.SetB(25, false);
	//std::cout << "Bc2 = " << get_bc2<t_real>(-1000.) << std::endl;

	std::cout.precision(8);
	std::cout << "Order: " << ORDER << std::endl;
	std::cout << "F_start = " << heli.F() << std::endl;
	bool ok = heli.minimise(ORDER, 0,1,0, 0,1,1);
	std::cout << "F_min = " << heli.F() << " (ok: " << std::boolalpha << ok << ")" << std::endl;

	std::cout << "\nFourier components:\n";
	for(const auto& fourier : heli.GetFourier())
	{
		std::cout
			<< "{ " << fourier[0].real() << ", " << fourier[0].imag() << " }, "
			<< "{ " << fourier[1].real() << ", " << fourier[1].imag() << " }, "
			<< "{ " << fourier[2].real() << ", " << fourier[2].imag() << " },"
			<< std::endl;
	}

	return 0;
}
