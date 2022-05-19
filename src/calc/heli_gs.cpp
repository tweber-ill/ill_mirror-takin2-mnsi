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
	Heli<t_real, t_cplx, DEF_HELI_ORDER> heli;

	t_real T_theo = -1000;
	t_real T_exp = 28.5;
	t_real B_theo = heli.GetBC2(false)/2.;
	t_real m_scale = 10.;

	//t_real T_theo = -4500;
	//t_real T_exp = 20.;
	//t_real B_theo = get_B_exp_from_theo(T_theo, T_exp, 0.171, !HELI_USE_HOC);
	//t_real m_scale = 15.;

	const auto j = t_cplx(0, 1);
	std::vector<ublas::vector<t_cplx>> fourier{
		tl2::make_vec<ublas::vector<t_cplx>>({0, 0, m_scale}),
		// helical order => Re{M} perp. Im{M}
		tl2::make_vec<ublas::vector<t_cplx>>({1.+j, 1.-j, 0}) / std::sqrt(2) * m_scale,
 	};

	heli.SetFourier(fourier);
	heli.SetT(T_theo, false);
	heli.SetB(B_theo, false);

	std::cout.precision(8);
	std::cout << "Order: " << ORDER << std::endl;
	std::cout << "T_theo = " << T_theo << ", T_exp = " << T_exp << std::endl;
	std::cout << "B_theo = " << B_theo << std::endl;
	std::cout << "F_start = " << heli.F() << std::endl;
	bool ok = heli.minimise(ORDER, 0,1,0, 0,1,1);
	std::cout << "F_min = " << heli.F() << " (ok: " << std::boolalpha << ok << ")" << std::endl;

	std::cout << "\nFourier components:\n";
	const auto& fouriers = heli.GetFourier();
	for(std::size_t pk_idx=0; pk_idx<fouriers.size(); ++pk_idx)
	{
		const auto& fourier = fouriers[pk_idx];
		int pk_num = abs_to_rel_idx((int)pk_idx, ORDER);

		std::cout
			<< "{ " << fourier[0].real() << ", " << fourier[0].imag() << " }, "
			<< "{ " << fourier[1].real() << ", " << fourier[1].imag() << " }, "
			<< "{ " << fourier[2].real() << ", " << fourier[2].imag() << " }, // " << pk_num
			<< std::endl;
	}

	return 0;
}
