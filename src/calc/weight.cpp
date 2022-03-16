/**
 * Calculates the weights/energies at a specific q
 * @author Tobias Weber <tweber@ill.fr>
 * @date sep-19
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/heli.h"
#include "core/skx.h"
#include "core/fp.h"

#include <fstream>
#include <memory>

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_vec = ublas::vector<t_real>;

#include "core/skx_default_gs.cxx"


void calc_disp(char dyntype,
	t_real Gx, t_real Gy, t_real Gz,
	t_real Bx, t_real By, t_real Bz,
	t_real Px, t_real Py, t_real Pz,
	t_real Qx, t_real Qy, t_real Qz,
	int iProj=1, t_real T=-1., t_real B=-1.)
{
	std::shared_ptr<MagDynamics<t_real, t_cplx>> dyn;

	if(dyntype == 's')
	{
		auto skx = std::make_shared<Skx<t_real, t_cplx, DEF_SKX_ORDER>>();

		std::vector<ublas::vector<t_cplx>> fourier_skx;
		fourier_skx.reserve(_skxgs_allcomps.size()/3);

		for(std::size_t comp=0; comp<_skxgs_allcomps.size(); comp+=3)
			fourier_skx.push_back(tl2::make_vec<ublas::vector<t_cplx>>({_skxgs_allcomps[comp], _skxgs_allcomps[comp+1], _skxgs_allcomps[comp+2]}));

		skx->SetFourier(fourier_skx);
		skx->SetProjNeutron(iProj!=0);
		skx->SetT(-1000.);
		skx->SetB(25.);
		skx->GenFullFourier();
		skx->SetFilterZeroWeight(1);
		skx->SetCoords(Bx,By,Bz, Px,Py,Pz);
		skx->SetG(Gx,Gy,Gz);

		dyn = skx;
	}
	else if(dyntype == 'h')
	{
		auto heli = std::make_shared<Heli<t_real, t_cplx, DEF_HELI_ORDER>>();

		heli->SetProjNeutron(iProj!=0);
		heli->SetT(T);
		heli->SetB(B);
		heli->SetFilterZeroWeight(1);
		heli->SetCoords(Bx,By,Bz);
		heli->SetG(Gx,Gy,Gz);

		dyn = heli;
	}
	else if(dyntype == 'f')
	{
		auto fp = std::make_shared<FP<t_real, t_cplx>>();

		fp->SetT(T);
		fp->SetB(B);
		fp->SetCoords(Bx,By,Bz);
		fp->SetG(Gx,Gy,Gz);

		dyn = fp;
	}
	else
	{
		std::cerr << "Unknown dynamics type selected." << std::endl;
		return;
	}


	t_real Erange = -1.;	// negative: disable range
	auto [Es, wsUnpol, wsSF1, wsSF2, wsNSF] = dyn->GetDisp(Qx, Qy, Qz, -Erange, Erange);

	std::cout << "# Magnetic phase: " << dyntype << "\n"
		<< std::setw(15) << std::left << "# No." << " "
		<< std::setw(15) << std::left << "E (meV)" << " "
		<< std::setw(15) << std::left << "w_total" << " "
		<< std::setw(15) << std::left << "w_SF1" << " "
		<< std::setw(15) << std::left << "w_SF2" << " "
		<< std::setw(15) << std::left << "w_NSF" << "\n";

	for(std::size_t i=0; i<Es.size(); ++i)
	{
		std::cout
			<< std::setw(15) << std::left << (i+1) << " "
			<< std::setw(15) << std::left /*<< std::scientific*/ << Es[i] << " "
			<< std::setw(15) << std::left /*<< std::scientific*/ << wsUnpol[i] << " "
			<< std::setw(15) << std::left /*<< std::scientific*/ << wsSF1[i] << " "
			<< std::setw(15) << std::left /*<< std::scientific*/ << wsSF2[i] << " "
			<< std::setw(15) << std::left /*<< std::scientific*/ << wsNSF[i] << "\n";
	}
}


int main()
{
	std::cout
		<< "--------------------------------------------------------------------------------\n"
		<< "\tDynamical structure factor calculation tool,\n\t\tT. Weber <tweber@ill.fr>, September 2019.\n"
		<< "--------------------------------------------------------------------------------\n\n";

	while(1)
	{
		t_real Gx = 1., Gy = 1., Gz = 0.;
		t_real Bx = 1., By = 1., Bz = 0.;
		t_real Qx = 1.068, Qy = 0.932, Qz = 0.;
		t_real Px = 1., Py = -1., Pz = 0.;
		t_real B = 0.17, T = 20.;
		int proj = 1;
		char dyntype = 'h';

		std::cout << "Helimagnon [h], skyrmion [s] or field-polarised [f] dynamics: ";
		std::cin >> dyntype; dyntype = std::tolower(dyntype);
		std::cout << "G = ";
		std::cin >> Gx >> Gy >> Gz;
		std::cout << "Q = ";
		std::cin >> Qx >> Qy >> Qz;
		std::cout << "B = ";
		std::cin >> Bx >> By >> Bz;
		if(dyntype == 'h' || dyntype == 'f')
		{
			std::cout << "|B| = ";
			std::cin >> B;
			std::cout << "T = ";
			std::cin >> T;
		}
		else if(dyntype == 's')
		{
			std::cout << "pinning = ";
			std::cin >> Px >> Py >> Pz;
		}
		std::cout << "Q projector [0/1]: ";
		std::cin >> proj;

		std::cout << "\n";
		if(dyntype == 'h' || dyntype == 'f')
		{
			std::cout << "# T = " << T << "\n";
			std::cout << "# |B| = " << B << "\n";
		}
		std::cout << "# B = (" << Bx << ", " << By << ", " << Bz << ")\n";
		std::cout << "# G = (" << Gx << ", " << Gy << ", " << Gz << ")\n";
		std::cout << "# Q = (" << Qx << ", " << Qy << ", " << Qz << ")\n";
		std::cout << "# q = (" << Qx-Gx << ", " << Qy-Gy << ", " << Qz-Gz << ")\n";
		std::cout << "# Q_proj = " << proj << "\n";

		calc_disp(dyntype, Gx,Gy,Gz, Bx,By,Bz, Px,Py,Pz, Qx,Qy,Qz, proj, T, B);
		std::cout << std::endl;
	}

	return 0;
}
