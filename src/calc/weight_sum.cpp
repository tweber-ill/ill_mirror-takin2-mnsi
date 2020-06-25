/**
 * Calculates the integrated weights/energies for the reseda experiment
 * @author Tobias Weber <tweber@ill.fr>
 * @date jun-20
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/heli.h"
#include "core/skx.h"

#include <fstream>
#include <future>
#include <pwd.h>

#include <boost/histogram.hpp>
namespace hist = boost::histogram;

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_vec = ublas::vector<t_real>;
const auto g_j = t_cplx(0,1);

#include "core/skx_default_gs.cxx"


#define E_BINS 200

void calc_disp(
	t_real Gx, t_real Gy, t_real Gz,
	t_real Bx, t_real By, t_real Bz,
	t_real Px, t_real Py, t_real Pz,
	t_real q,
	int iProj=1)
{
	Skx<t_real, t_cplx, DEF_SKX_ORDER> skx;

	std::vector<ublas::vector<t_cplx>> fourier_skx;
	fourier_skx.reserve(_skxgs_allcomps.size()/3);

	for(std::size_t comp=0; comp<_skxgs_allcomps.size(); comp+=3)
		fourier_skx.push_back(tl2::make_vec<ublas::vector<t_cplx>>({_skxgs_allcomps[comp], _skxgs_allcomps[comp+1], _skxgs_allcomps[comp+2]}));

	skx.SetFourier(fourier_skx);
	skx.SetProjNeutron(iProj!=0);
	skx.SetT(-1000.);
	skx.SetB(25.);	// BC2 = 40.3425
	skx.GenFullFourier();
 	skx.SetFilterZeroWeight(1);
	skx.SetWeightEps(1e-6);

	skx.SetCoords(Bx,By,Bz, Px,Py,Pz);
	skx.SetG(Gx, Gy, Gz);
	t_real Erange = 0.1;

	t_real angle_begin = -135/180.*tl2::pi<t_real>;
	t_real angle_end = 135/180.*tl2::pi<t_real>;
	t_real angle_delta = 2*tl2::pi<t_real>/100.;

	auto histWeights = hist::make_histogram(hist::axis::regular<t_real>(E_BINS, -Erange, Erange, "E"));

	for(t_real angle=angle_begin; angle<angle_end; angle+=angle_delta)
	{
		t_real Qx = q * cos(angle) + Gx;
		t_real Qy = q * sin(angle) + Gy;
		t_real Qz = Gz;

		std::cout << "# angle: " << angle/tl2::pi<t_real>*180. << ", Q = (" << Qx << ", " << Qy << ", " << Qz << ")\n";

		auto [Es, wsUnpol, wsSF1, wsSF2, wsNSF] = skx.GetDisp(Qx, Qy, Qz, -Erange, Erange);
		for(std::size_t i=0; i<Es.size(); ++i)
			histWeights(Es[i], hist::weight(wsNSF[i]*0.5 + wsSF1[i]));
	}


	std::cout << std::left << std::setw(15) << "# E (meV)" << " " << std::left << std::setw(15) << "weight" << "\n";

	for(const auto& val : boost::histogram::indexed(histWeights))
	{
		t_real E = val.bin().lower() + 0.5*(val.bin().upper() - val.bin().lower());
		t_real w = *val;

		std::cout << std::left << std::setw(15) << E << " " << std::left << std::setw(15) << w << "\n";
	}
}


int main()
{
	std::cout.precision(5);

	// kh = 0.02829 rlu
	t_real Gx = 1., Gy = 0., Gz = 0.;
	t_real Bx = 0., By = 0., Bz = 1.;
	t_real q = 0.0123;
	int proj = 1;

	// setup used at reseda
	calc_disp(Gx,Gy,Gz, Bx,By,Bz, 1,0,0, q, proj);

	return 0;
}
