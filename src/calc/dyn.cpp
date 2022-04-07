/**
 * Calculates dispersion curves
 * @author Tobias Weber <tweber@ill.fr>
 * @date aug-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/heli.h"
#include "core/skx.h"
#include "core/fp.h"

#include <fstream>
#include <future>
#include <memory>

#ifndef __MINGW32__
	#include <pwd.h>
#endif

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_vec = ublas::vector<t_real>;
using t_vec_cplx = ublas::vector<t_cplx>;

#include "core/heli_default_gs.cxx"
#include "core/skx_default_gs.cxx"


static void print_groundstate(const std::vector<t_vec_cplx>& gs)
{
	std::cout << "Ground state:\n";
	for(const auto& fourier : gs)
	{
		std::cout
			<< "{ " << fourier[0].real() << ", " << fourier[0].imag() << " }, "
			<< "{ " << fourier[1].real() << ", " << fourier[1].imag() << " }, "
			<< "{ " << fourier[2].real() << ", " << fourier[2].imag() << " },"
			<< std::endl;
	}
}


static void calc_disp(char dyntype,
	t_real Gx, t_real Gy, t_real Gz,
	t_real Bx, t_real By, t_real Bz,
	t_real Px, t_real Py, t_real Pz,
	t_real qperpdir_x, t_real qperpdir_y, t_real qperpdir_z, t_real qperp,
	const char* pcFile, bool bSwapQParaQPerp=0,
	t_real T=-1., t_real B=-1.,
	t_real qrange = 0.125, t_real delta = 0.001)
{
	tl2::Stopwatch<t_real> timer;
	timer.start();

	constexpr auto imag = t_cplx(0, 1);
	t_vec G = tl2::make_vec<t_vec>({ Gx, Gy, Gz });
	t_vec Pdir = tl2::make_vec<t_vec>({ Px, Py, Pz });
	t_vec Bdir = tl2::make_vec<t_vec>({ Bx, By, Bz });
	t_vec qparadir = Bdir / tl2::veclen(Bdir);
	t_vec qperpdir = tl2::make_vec<t_vec>({ qperpdir_x, qperpdir_y, qperpdir_z });
	qperpdir /= tl2::veclen(qperpdir);

	std::shared_ptr<MagDynamics<t_real, t_cplx>> dyn;

	if(dyntype == 's')
	{
		std::cout << "Calculating skyrmion dispersion." << std::endl;
		auto skx = std::make_shared<Skx<t_real, t_cplx, DEF_SKX_ORDER>>();
		skx->SetFourier(_get_skx_gs<t_vec_cplx>());
		print_groundstate(skx->GetFourier());
		dyn = skx;
	}
	else if(dyntype == 'h')
	{
		std::cout << "Calculating helical dispersion." << std::endl;
		auto heli = std::make_shared<Heli<t_real, t_cplx, DEF_HELI_ORDER>>();
		heli->SetFourier(_get_heli_gs<t_vec_cplx>());
		print_groundstate(heli->GetFourier());
		dyn = heli;

	}
	else if(dyntype == 'f')
	{
		std::cout << "Calculating field-polarised dispersion." << std::endl;
		dyn = std::make_shared<FP<t_real, t_cplx>>();
	}
	else
	{
		std::cerr << "Unknown dynamics type selected." << std::endl;
		return;
	}


	dyn->SetCoords(Bdir[0],Bdir[1],Bdir[2], Pdir[0],Pdir[1],Pdir[2]);
	dyn->SetT(-1000., false);
	dyn->SetB(25., false);	// BC2 = 45.028487
	dyn->SetT(T, true);
	dyn->SetB(B, true);
	dyn->SetG(G[0], G[1], G[2]);


	t_real F = 0.;
	auto *magsys = dynamic_cast<HasF<t_real>*>(dyn.get());
	if(magsys)
	{
		F = magsys->F();
		std::cout << "Ground state F = " << F << "." << std::endl;
	}


	auto calc_spectrum = [dyntype, &dyn, &G, T, &qparadir, &qperpdir, &qperp, bSwapQParaQPerp]
		(int thid, t_real qstart, t_real qend, t_real qdelta) -> auto
	{
		std::vector<t_real> allh, allk, alll, allqpara_kh, allqperp_kh, allqpara_rlu, allqperp_rlu;
		std::vector<std::vector<t_real>> allEs, allWsUnpol, allWsSF1, allWsSF2, allWsNSF;
		auto thisdyn = dyn->copyCastDyn();

		for(t_real _q=qstart; _q<qend; _q+=qdelta)
		{
			t_real qpara = _q;
			t_vec Q = G + qpara*qparadir + qperp*qperpdir;
			if(bSwapQParaQPerp)
				Q = G + qpara*qperpdir + qperp*qparadir;	// swap qpara and qperp

			std::cout << "thread " << thid << " (" << Q[0] << " " << Q[1] << " " << Q[2] << ") ... ";
			std::cout.flush();

			auto [Es, wsUnpol, wsSF1, wsSF2, wsNSF] = thisdyn->GetDisp(Q[0], Q[1], Q[2]);
			allEs.emplace_back(Es);
			allWsUnpol.emplace_back(wsUnpol);
			allWsSF1.emplace_back(wsSF1);
			allWsSF2.emplace_back(wsSF2);
			allWsNSF.emplace_back(wsNSF);

			allh.push_back(Q[0]);
			allk.push_back(Q[1]);
			alll.push_back(Q[2]);

			if(dyntype == 's')
			{
				allqpara_kh.push_back(qpara / g_kh_rlu_29K<t_real>);
				allqperp_kh.push_back(qperp / g_kh_rlu_29K<t_real>);
			}
			else
			{
				allqpara_kh.push_back(qpara / g_kh_rlu<t_real>(T));
				allqperp_kh.push_back(qperp / g_kh_rlu<t_real>(T));
			}

			allqpara_rlu.push_back(qpara);
			allqperp_rlu.push_back(qperp);

			std::cout << "done." << std::endl;
		}

		return std::make_tuple(allh, allk, alll, allEs,
			allWsUnpol, allWsSF1, allWsSF2, allWsNSF,
			allqpara_kh, allqperp_kh, allqpara_rlu, allqperp_rlu);
	};


	t_real qstep = qrange / 2.;

	auto fut0 = std::async(std::launch::async, calc_spectrum, 0, -qrange, -qrange+1.*qstep, delta);
	auto fut1 = std::async(std::launch::async, calc_spectrum, 1, -qrange+1.*qstep, -qrange+2.*qstep, delta);
	auto fut2 = std::async(std::launch::async, calc_spectrum, 2, -qrange+2.*qstep, -qrange+3.*qstep, delta);
	auto fut3 = std::async(std::launch::async, calc_spectrum, 3, -qrange+3.*qstep, -qrange+4.*qstep, delta);
	auto val0 = fut0.get();
	auto val1 = fut1.get();
	auto val2 = fut2.get();
	auto val3 = fut3.get();

	insert_vals(val0, val1, std::make_index_sequence<std::tuple_size<decltype(val1)>::value>());
	insert_vals(val0, val2, std::make_index_sequence<std::tuple_size<decltype(val2)>::value>());
	insert_vals(val0, val3, std::make_index_sequence<std::tuple_size<decltype(val3)>::value>());


	timer.stop();
	std::ofstream ofstr(pcFile);
	ofstr.precision(8);

	ofstr << "#\n";
	ofstr << "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n";
#ifndef __MINGW32__
	ofstr << "# User: " << getpwuid(geteuid())->pw_name << "\n";
#endif
	ofstr << "# Calculation time: " << timer.GetDur() << " s.\n";
	ofstr << "# Skx order: " << DEF_SKX_ORDER << "\n";
	if(dyntype == 'h')
	{
		ofstr << "# Heli order: " << DEF_HELI_ORDER << "\n";
		ofstr << "# kh_A = " << g_kh_A<t_real>(T) << "\n";
		ofstr << "# kh_rlu = " << g_kh_rlu<t_real>(T) << "\n";
	}
	ofstr << "# qparadir = " << qparadir << "\n";
	ofstr << "# qperpdir = " << qperpdir << "\n";
	ofstr << "# qperp = " << qperp << "\n";
	ofstr << "# G = " << G << "\n";
	ofstr << "# Bdir = " << Bdir << "\n";
	ofstr << "# Pdir = " << Pdir << "\n";
	ofstr << "# Bc2 = " << dyn->GetBC2(false) << "\n";
	ofstr << "# Bc2_exp = " << dyn->GetBC2(true) << "\n";
	ofstr << "# F = " << F << "\n";
	if(bSwapQParaQPerp)
		ofstr << "# WARNING: In the following, q_para_* and q_perp_* are swapped!\n";

	// save commutator of projectors to determine channel mixing
	/*ofstr << "# [Nrpoj, Polproj_sf1] = " << tl2::commutator<t_mat_cplx>(fp.GetNeutronProjOp(), get_chiralpol<t_mat_cplx>(1)) << "\n";
	ofstr << "# [Nrpoj, Polproj_sf2] = " << tl2::commutator<t_mat_cplx>(fp.GetNeutronProjOp(), get_chiralpol<t_mat_cplx>(2)) << "\n";
	ofstr << "# [Nrpoj, Polproj_nsf] = " << tl2::commutator<t_mat_cplx>(fp.GetNeutronProjOp(), get_chiralpol<t_mat_cplx>(3)) << "\n";*/

	ofstr << "#\n";

	ofstr
		<< "#" << std::setw(15) << "h" << " "
		<< std::setw(16) << "k" << " "
		<< std::setw(16) << "l" << " "
		<< std::setw(16) << "E" << " "
		<< std::setw(16) << "w_unpol" << " "
		<< std::setw(16) << "w_sf1" << " "
		<< std::setw(16) << "w_sf2" << " "
		<< std::setw(16) << "w_nsf" << " "
		<< std::setw(16) << "q_para_kh" << " "
		<< std::setw(16) << "q_perp_kh" << " "
		<< std::setw(16) << "q_para_rlu" << " "
		<< std::setw(16) << "q_perp_rlu\n";

	for(std::size_t i=0; i<std::get<0>(val0).size(); ++i)
	{
		for(std::size_t j=0; j<std::get<3>(val0)[i].size(); ++j)
		{
			ofstr << std::setw(16) << std::get<0>(val0)[i] << " "       // h
				<< std::setw(16) << std::get<1>(val0)[i] << " "     // k
				<< std::setw(16) << std::get<2>(val0)[i] << " "     // l
				<< std::setw(16) << std::get<3>(val0)[i][j] << " "  // E
				<< std::setw(16) << std::get<4>(val0)[i][j] << " "  // w_unpol
				<< std::setw(16) << std::get<5>(val0)[i][j] << " "  // w_sf1
				<< std::setw(16) << std::get<6>(val0)[i][j] << " "  // w_sf2
				<< std::setw(16) << std::get<7>(val0)[i][j] << " "  // w_nsf
				<< std::setw(16) << std::get<8>(val0)[i] << " "     // q_para_kh
				<< std::setw(16) << std::get<9>(val0)[i] << " "     // q_perp_kh
				<< std::setw(16) << std::get<10>(val0)[i] << " "    // q_para_rlu
				<< std::setw(16) << std::get<11>(val0)[i] << "\n";  // q_perp_rlu
		}
	}


	timer.stop();
	std::cout << "Calculation took " << timer.GetDur() << " s." << std::endl;
	std::cout << "Wrote \"" << pcFile << "\"" << std::endl;
}


int main()
{
	std::cout
		<< "--------------------------------------------------------------------------------\n"
		<< "\tDispersion calculation tool,\n\t\tT. Weber <tweber@ill.fr>, August 2018.\n"
		<< "--------------------------------------------------------------------------------\n\n";

	char dyntype = 's';
	t_real Gx = 1., Gy = 1., Gz = 0.;
	t_real Bx = 1., By = 1., Bz = 0.;
	t_real Px = -1., Py = 1., Pz = 0.;
	t_real qperpx = -1., qperpy = 1., qperpz = 0.;
	t_real qperp = 0.;
	t_real B = 0.17, T = 28.5;
	t_real qrange = 0.125;
	t_real qdelta = 0.001;
	int alongqpara = 0;

	std::cout << "Helimagnon [h], skyrmion [s] or field-polarised [f] dynamics: ";
	std::cin >> dyntype; dyntype = std::tolower(dyntype);
	std::cout << "G = ";
	std::cin >> Gx >> Gy >> Gz;
	std::cout << "q_range = ";
	std::cin >> qrange;
	std::cout << "q_delta = ";
	std::cin >> qdelta;
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
	std::cout << "Query dispersion along q_para || B? [1/0]: ";
	std::cin >> alongqpara;
	std::cout << "q_perp = ";
	std::cin >> qperpx >> qperpy >> qperpz;
	std::cout << "|q_perp| = ";
	std::cin >> qperp;

	calc_disp(dyntype, Gx,Gy,Gz, Bx,By,Bz, Px,Py,Pz,
		qperpx,qperpy,qperpz, qperp,
		"dyn.dat", alongqpara==0,
		T, B,
		qrange, qdelta);

	return 0;
}
