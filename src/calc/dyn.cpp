/**
 * Calculates dispersion curves
 * @author Tobias Weber <tweber@ill.fr>
 * @date aug-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/heli.h"
#include "core/skx.h"
#include "core/fp.h"
#include "core/load_gs.h"

#include <fstream>
#include <thread>
#include <future>
#include <memory>

#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>
#include <boost/asio.hpp>

namespace asio = boost::asio;
namespace opts = boost::program_options;

#ifndef __MINGW32__
	#include <pwd.h>
#endif

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_vec = ublas::vector<t_real>;
using t_vec_cplx = ublas::vector<t_cplx>;
using t_mat = ublas::matrix<t_real>;

#include "core/heli_default_gs.cxx"
#include "core/skx_default_gs.cxx"


/**
 * print ground state magnetisation
 */
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


static void show_progress(t_real done, bool init = false)
{
	static int last_percentage = -1;
	int new_percentage = static_cast<int>(done * 100.);
	if(init)
		last_percentage = -1;

	if(new_percentage != last_percentage)
	{
		std::cout << "\rRunning, " << new_percentage << "% done.        ";
		std::cout.flush();

		last_percentage = new_percentage;
	}
}


/**
 * create the skyrmion, conical, or field-polarised dynamics module
 */
static std::shared_ptr<MagDynamics<t_real, t_cplx>> create_dyn(
	char dyntype, bool explicit_calc = true, const std::string& gs_file = "")
{
	std::shared_ptr<MagDynamics<t_real, t_cplx>> dyn;

	if(dyntype == 's')
	{
		//std::cout << "Calculating skyrmion dispersion." << std::endl;

		// get default ground state
		auto [skxgs_T, skxgs_B, skxgs] = _get_skx_gs<t_vec_cplx>();

		// optionally load a given ground state file
		if(gs_file != "")
		{
			bool ok = false;
			std::tie(ok, skxgs_T, skxgs_B, skxgs) =
				load_gs<std::decay_t<decltype(skxgs)>>(gs_file, 's');
			if(!ok)
			{
				std::cerr << "Error: Could not load skyrmion ground state." << std::endl;
				return nullptr;
			}
		}

		auto skx = std::make_shared<Skx<t_real, t_cplx, DEF_SKX_ORDER>>();
		skx->SetFourier(skxgs);
		skx->SetT(skxgs_T, false);
		skx->SetB(/*skx->GetBC2(false)/2.*/ skxgs_B, false);

		print_groundstate(skx->GetFourier());
		dyn = skx;
	}
	else if(dyntype == 'h')
	{
		//std::cout << "Calculating helical dispersion." << std::endl;

		// get default ground state
		auto [heligs_T, heligs_B, heligs] = _get_heli_gs<t_vec_cplx>();

		// optionally load a given ground state file
		if(gs_file != "")
		{
			bool ok = false;
			std::tie(ok, heligs_T, heligs_B, heligs) =
				load_gs<std::decay_t<decltype(heligs)>>(gs_file, 'h');
			if(!ok)
			{
				std::cerr << "Error: Could not load conical ground state." << std::endl;
				return nullptr;
			}
		}

		auto heli = std::make_shared<Heli<t_real, t_cplx, DEF_HELI_ORDER>>();
		heli->SetExplicitCalc(explicit_calc);
		heli->SetFourier(heligs);
		heli->SetT(heligs_T, false);
		heli->SetB(/*heli->GetBC2(false)/2.*/ heligs_B, false);

		print_groundstate(heli->GetFourier());
		dyn = heli;
	}
	else if(dyntype == 'f')
	{
		//std::cout << "Calculating field-polarised dispersion." << std::endl;
		dyn = std::make_shared<FP<t_real, t_cplx>>();
	}
	else
	{
		std::cerr << "Error: Unknown dynamics type selected." << std::endl;
		return nullptr;
	}

	return dyn;
}


/**
 * calculate the dispersion parallel and perpendicular to the field
 */
static void calc_disp_para_perp(char dyntype, bool do_proj,
	t_real Gx, t_real Gy, t_real Gz,
	t_real Bx, t_real By, t_real Bz,
	t_real Px, t_real Py, t_real Pz,
	t_real qperpdir_x, t_real qperpdir_y, t_real qperpdir_z,
	t_real qperp, t_real qperp2,
	const std::string& outfile, bool bSwapQParaQPerp=0,
	t_real T=-1., t_real B=-1.,
	t_real qrange = 0.125, t_real delta = 0.001,
	bool explicit_calc = true,
	const std::string& gs_file = "")
{
	tl2::Stopwatch<t_real> timer;
	timer.start();

	std::shared_ptr<MagDynamics<t_real, t_cplx>> dyn =
		create_dyn(dyntype, explicit_calc, gs_file);
	if(!dyn)
		return;

	t_vec G = tl2::make_vec<t_vec>({ Gx, Gy, Gz });
	t_vec Pdir = tl2::make_vec<t_vec>({ Px, Py, Pz });
	t_vec Bdir = tl2::make_vec<t_vec>({ Bx, By, Bz });

	// direction parallel to the field
	t_vec qparadir = Bdir / tl2::veclen(Bdir);

	// direction perpendicular to the field
	t_vec qperpdir = tl2::make_vec<t_vec>({ qperpdir_x, qperpdir_y, qperpdir_z });
	qperpdir /= tl2::veclen(qperpdir);

	t_vec qperpdir2 = tl2::cross_3(qparadir, qperpdir);
	qperpdir2 /= tl2::veclen(qperpdir2);

	dyn->SetCoords(Bdir[0],Bdir[1],Bdir[2], Pdir[0],Pdir[1],Pdir[2]);
	dyn->SetT(T, true);
	dyn->SetB(B, true);
	dyn->SetProjNeutron(do_proj);
	dyn->SetG(G[0], G[1], G[2]);


	t_real F = 0.;
	auto *magsys = dynamic_cast<HasF<t_real>*>(dyn.get());
	if(magsys)
	{
		F = magsys->F();
		std::cout << "Ground state F = " << F << "." << std::endl;
	}

	std::cout << "Bc2_theo = " << dyn->GetBC2(false)
		<< ", Bc2_exp = " << dyn->GetBC2(true)
		<< "." << std::endl;


	auto calc_spectrum = [dyntype, &dyn, &G, T, &qparadir,
		&qperpdir, &qperp, &qperpdir2, &qperp2, bSwapQParaQPerp]
		(int thid, t_real qstart, t_real qend, t_real qdelta) -> auto
	{
		std::vector<t_real> allh, allk, alll, allqpara_kh, allqperp_kh, allqpara_rlu, allqperp_rlu;
		std::vector<std::vector<t_real>> allEs, allWsUnpol, allWsSF1, allWsSF2, allWsNSF;
		auto thisdyn = dyn->copyCastDyn();

		for(t_real _q=qstart; _q<qend; _q+=qdelta)
		{
			t_real qpara = _q;
			t_vec Q = G + qpara*qparadir + qperp*qperpdir + qperp2*qperpdir2;
			if(bSwapQParaQPerp)
				Q = G + qpara*qperpdir + qperp*qparadir + qperp2*qperpdir2;	// swap qpara and qperp

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


	// write results
	std::ofstream ofstr(outfile);
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
			ofstr << std::setw(16) << std::get<0>(val0)[i] << " "       // 1: h
				<< std::setw(16) << std::get<1>(val0)[i] << " "     // 2: k
				<< std::setw(16) << std::get<2>(val0)[i] << " "     // 3: l
				<< std::setw(16) << std::get<3>(val0)[i][j] << " "  // 4: E
				<< std::setw(16) << std::get<4>(val0)[i][j] << " "  // 5: w_unpol
				<< std::setw(16) << std::get<5>(val0)[i][j] << " "  // 6: w_sf1
				<< std::setw(16) << std::get<6>(val0)[i][j] << " "  // 7: w_sf2
				<< std::setw(16) << std::get<7>(val0)[i][j] << " "  // 8: w_nsf
				<< std::setw(16) << std::get<8>(val0)[i] << " "     // 9: q_para_kh
				<< std::setw(16) << std::get<9>(val0)[i] << " "     // 10: q_perp_kh
				<< std::setw(16) << std::get<10>(val0)[i] << " "    // 11: q_para_rlu
				<< std::setw(16) << std::get<11>(val0)[i] << "\n";  // 12: q_perp_rlu
		}
	}

	std::cout << "Calculation took " << timer.GetDur() << " s." << std::endl;
	std::cout << "Wrote \"" << outfile << "\"." << std::endl;
}


/**
 * calculate the dispersion along arbitrary paths of momentum transfer
 */
static void calc_disp_path(char dyntype, bool do_proj,
	t_real Gx, t_real Gy, t_real Gz,
	t_real Bx, t_real By, t_real Bz,
	t_real Px, t_real Py, t_real Pz,
	t_real qh_start, t_real qk_start, t_real ql_start,
	t_real qh_end, t_real qk_end, t_real ql_end,
	t_real Rx, t_real Ry, t_real Rz, t_real Ralpha,
	std::size_t num_points, const std::string& outfile,
	t_real T=-1., t_real B=-1., bool explicit_calc = true,
	const std::string& gs_file = "",
	unsigned int num_threads = 4)
{
	tl2::Stopwatch<t_real> timer;
	timer.start();

	std::shared_ptr<MagDynamics<t_real, t_cplx>> dyn =
		create_dyn(dyntype, explicit_calc, gs_file);
	if(!dyn)
		return;

	t_vec G = tl2::make_vec<t_vec>({ Gx, Gy, Gz });
	t_vec Pdir = tl2::make_vec<t_vec>({ Px, Py, Pz });
	t_vec Bdir = tl2::make_vec<t_vec>({ Bx, By, Bz });

	t_vec qstart = tl2::make_vec<t_vec>({ qh_start, qk_start, ql_start });
	t_vec qend = tl2::make_vec<t_vec>({ qh_end, qk_end, ql_end });

	t_vec Rdir = tl2::make_vec<t_vec>({ Rx, Ry, Rz });
	t_mat rot = tl2::rotation_matrix<t_mat, t_vec>(Rdir, Ralpha);
	qstart = tl2::prod_mv(rot, qstart);
	qend = tl2::prod_mv(rot, qend);

	dyn->SetCoords(Bdir[0],Bdir[1],Bdir[2], Pdir[0],Pdir[1],Pdir[2]);
	dyn->SetT(T, true);
	dyn->SetB(B, true);
	dyn->SetProjNeutron(do_proj);
	dyn->SetG(G[0], G[1], G[2]);

	t_real F = 0.;
	auto *magsys = dynamic_cast<HasF<t_real>*>(dyn.get());
	if(magsys)
	{
		F = magsys->F();
		std::cout << "Ground state F = " << F << "." << std::endl;
	}

	std::cout << "Bc2_theo = " << dyn->GetBC2(false)
		<< ", Bc2_exp = " << dyn->GetBC2(true)
		<< "." << std::endl;

	// thread pool
	asio::thread_pool pool(num_threads);

	// calculation task
	using t_task = std::packaged_task<void()>;
	std::vector<std::shared_ptr<t_task>> tasks;
	tasks.reserve(num_points);

	std::vector<t_real> allh(num_points), allk(num_points), alll(num_points);
	std::vector<std::vector<t_real>> allEs(num_points),
		allWsUnpol(num_points), allWsSF1(num_points),
		allWsSF2(num_points), allWsNSF(num_points);
	std::vector<t_real> allqs(num_points);

	for(std::size_t pt_idx=0; pt_idx<num_points; ++pt_idx)
	{
		t_vec q = qstart + t_real(pt_idx)/t_real(num_points-1) * (qend-qstart);
		t_vec Q = G + q;

		t_real qlen = tl2::veclen(q);
		t_real qsign = tl2::inner<t_vec>(qend-qstart, q) < 0. ? -1. : 1.;

		auto calc_spectrum = [&dyn, pt_idx, Q, qsign, qlen,
			&allh, &allk, &alll, &allEs, &allWsUnpol,
			&allWsSF1, &allWsSF2, &allWsNSF,
			&allqs]()
		{
			//auto thisdyn = dyn->copyCastDyn();
			auto& thisdyn = dyn;
			auto [Es, wsUnpol, wsSF1, wsSF2, wsNSF] = thisdyn->GetDisp(Q[0], Q[1], Q[2]);

			allh[pt_idx] = Q[0]; allk[pt_idx] = Q[1]; alll[pt_idx] = Q[2];
			allEs[pt_idx] = Es;
			allWsUnpol[pt_idx] = wsUnpol;
			allWsSF1[pt_idx] = wsSF1;
			allWsSF2[pt_idx] = wsSF2;
			allWsNSF[pt_idx] = wsNSF;
			allqs[pt_idx] = qsign * qlen;
		};

		auto taskptr = std::make_shared<t_task>(calc_spectrum);
		tasks.push_back(taskptr);
		asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	show_progress(0, true);
	for(std::size_t taskidx=0; taskidx<tasks.size(); ++taskidx)
	{
		tasks[taskidx]->get_future().get();
		t_real done = t_real(taskidx+1) / t_real(tasks.size());
		show_progress(done);
	}

	pool.join();
	timer.stop();


	// write results
	std::ofstream ofstr(outfile);
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
	ofstr << "# G = " << G << "\n";
	ofstr << "# Bdir = " << Bdir << "\n";
	ofstr << "# Pdir = " << Pdir << "\n";
	ofstr << "# Bc2 = " << dyn->GetBC2(false) << "\n";
	ofstr << "# Bc2_exp = " << dyn->GetBC2(true) << "\n";
	ofstr << "# F = " << F << "\n";
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
		<< std::setw(16) << "q\n";

	for(std::size_t i=0; i<num_points; ++i)
	{
		for(std::size_t j=0; j<allEs[i].size(); ++j)
		{
			ofstr
				<< std::setw(16) << allh[i] << " "           // 1: h
				<< std::setw(16) << allk[i] << " "           // 2: k
				<< std::setw(16) << alll[i] << " "           // 3: l
				<< std::setw(16) << allEs[i][j] << " "       // 4: E
				<< std::setw(16) << allWsUnpol[i][j] << " "  // 5: w_unpol
				<< std::setw(16) << allWsSF1[i][j] << " "    // 6: w_sf1
				<< std::setw(16) << allWsSF2[i][j] << " "    // 7: w_sf2
				<< std::setw(16) << allWsNSF[i][j] << " "    // 8: w_nsf
				<< std::setw(16) << allqs[i] << "\n";        // 9: q
		}
	}

	std::cout << "\nCalculation took " << timer.GetDur() << " s." << std::endl;
	std::cout << "Wrote \"" << outfile << "\"." << std::endl;
}


int main(int argc, char **argv)
{
	std::cout
		<< "--------------------------------------------------------------------------------\n"
		<< "\tDispersion calculation tool,\n\t\tT. Weber <tweber@ill.fr>, August 2018.\n"
		<< "--------------------------------------------------------------------------------\n\n";

	// arguments for simple Q path calculation
	char dyntype = 's';
	t_real Gx = 1., Gy = 1., Gz = 0.;
	t_real Bx = 1., By = 1., Bz = 0.;
	t_real Px = -1., Py = 1., Pz = 0.;
	t_real qperpx = -1., qperpy = 1., qperpz = 0.;
	t_real qperp = 0., qperp2 = 0.;
	t_real B = 0.17, T = 28.5;
	t_real qrange = 0.125;
	t_real qdelta = 0.001;
	bool alongqpara = false;
	bool do_proj = true;
	bool explicit_calc = true;
	bool use_para_perp_calc = true;
	std::string outfile = "dyn.dat";
	std::string gs_file = "";

	// arguments for arbitrary Q path calculation
	t_real qh_start = 0., qk_start = 0., ql_start = 0.;
	t_real qh_end = 0.1, qk_end = 0., ql_end = 0.;
	t_real Rx = 0., Ry = 0., Rz = 1., Ralpha = 0.;
	std::size_t num_points = 256;
	unsigned int num_threads =
		std::max<unsigned int>(1, std::thread::hardware_concurrency()/2);


	if(argc <= 1)
	{
		std::cout << "No arguments given, running interactively. "
			<< "Use \"" << argv[0] << " --help\" to show all options.\n"
			<< std::endl;

		std::cout << "Helimagnon [h], skyrmion [s] or field-polarised [f] dynamics: ";
		std::cin >> dyntype;
		dyntype = std::tolower(dyntype);

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

			if(dyntype == 'h')
			{
				std::cout << "Explicit calculation? [1/0]: ";
				std::cin >> explicit_calc;
			}
		}
		else if(dyntype == 's')
		{
			std::cout << "pinning = ";
			std::cin >> Px >> Py >> Pz;
		}
		std::cout << "Calculate S_perp? [1/0]: ";
		std::cin >> do_proj;
		std::cout << "Query dispersion along q_para || B? [1/0]: ";
		std::cin >> alongqpara;
		std::cout << "q_perp = ";
		std::cin >> qperpx >> qperpy >> qperpz;
		std::cout << "|q_perp| = ";
		std::cin >> qperp;
		std::cout << "|q_perp_2| = ";
		std::cin >> qperp2;
	}
	else
	{
		try
		{
			opts::basic_command_line_parser<char> clparser(argc, argv);
			opts::options_description args("program arguments");

			bool show_help = false;
			args.add(boost::make_shared<opts::option_description>(
				"help", opts::bool_switch(&show_help),
				"show program usage"));

			args.add(boost::make_shared<opts::option_description>(
				"use_para_perp_calc", opts::value<decltype(use_para_perp_calc)>(&use_para_perp_calc),
				"use simple Q path calculation"));

			args.add(boost::make_shared<opts::option_description>(
				"dyntype", opts::value<decltype(dyntype)>(&dyntype),
				"dispersion type [s/h/f]"));

			args.add(boost::make_shared<opts::option_description>(
				"explicit_calc", opts::value<decltype(explicit_calc)>(&explicit_calc),
				"use explicit calculation"));

			args.add(boost::make_shared<opts::option_description>(
				"along_qpara", opts::value<decltype(alongqpara)>(&alongqpara),
				"calculate dispersion along q_parallel"));

			args.add(boost::make_shared<opts::option_description>(
				"do_proj", opts::value<decltype(do_proj)>(&do_proj),
				"calculate orthogonal projection, S_perp"));

			args.add(boost::make_shared<opts::option_description>(
				"outfile", opts::value<decltype(outfile)>(&outfile),
				"output file name"));

			args.add(boost::make_shared<opts::option_description>(
				"gsfile", opts::value<decltype(gs_file)>(&gs_file),
				"ground state file name"));

			args.add(boost::make_shared<opts::option_description>(
				"Gx", opts::value<decltype(Gx)>(&Gx),
				"lattice vector x component"));
			args.add(boost::make_shared<opts::option_description>(
				"Gy", opts::value<decltype(Gy)>(&Gy),
				"lattice vector y component"));
			args.add(boost::make_shared<opts::option_description>(
				"Gz", opts::value<decltype(Gz)>(&Gz),
				"lattice vector z component"));

			args.add(boost::make_shared<opts::option_description>(
				"Bx", opts::value<decltype(Bx)>(&Bx),
				"magnetic field vector x component"));
			args.add(boost::make_shared<opts::option_description>(
				"By", opts::value<decltype(By)>(&By),
				"magnetic field vector y component"));
			args.add(boost::make_shared<opts::option_description>(
				"Bz", opts::value<decltype(Bz)>(&Bz),
				"magnetic field vector z component"));

			args.add(boost::make_shared<opts::option_description>(
				"Px", opts::value<decltype(Px)>(&Px),
				"pinning vector x component"));
			args.add(boost::make_shared<opts::option_description>(
				"Py", opts::value<decltype(Py)>(&Py),
				"pinning vector y component"));
			args.add(boost::make_shared<opts::option_description>(
				"Pz", opts::value<decltype(Pz)>(&Pz),
				"pinning vector z component"));

			args.add(boost::make_shared<opts::option_description>(
				"qperpx", opts::value<decltype(qperpx)>(&qperpx),
				"perpendicular q direction x component"));
			args.add(boost::make_shared<opts::option_description>(
				"qperpy", opts::value<decltype(qperpy)>(&qperpy),
				"perpendicular q direction y component"));
			args.add(boost::make_shared<opts::option_description>(
				"qperpz", opts::value<decltype(qperpz)>(&qperpz),
				"perpendicular q direction z component"));
			args.add(boost::make_shared<opts::option_description>(
				"qperp", opts::value<decltype(qperp)>(&qperp),
				"perpendicular q magnitude"));
			args.add(boost::make_shared<opts::option_description>(
				"qperp2", opts::value<decltype(qperp2)>(&qperp2),
				"second perpendicular q magnitude"));

			args.add(boost::make_shared<opts::option_description>(
				"qrange", opts::value<decltype(qrange)>(&qrange),
				"parallel q range"));
			args.add(boost::make_shared<opts::option_description>(
				"qdelta", opts::value<decltype(qdelta)>(&qdelta),
				"parallel q delta"));

			args.add(boost::make_shared<opts::option_description>(
				"T", opts::value<decltype(T)>(&T),
				"temperature"));
			args.add(boost::make_shared<opts::option_description>(
				"B", opts::value<decltype(B)>(&B),
				"magnetic field magnitude"));


			// arguments for arbitrary Q path calculation
			args.add(boost::make_shared<opts::option_description>(
				"num_points", opts::value<decltype(num_points)>(&num_points),
				"number of points in Q path"));

			args.add(boost::make_shared<opts::option_description>(
				"qh_start", opts::value<decltype(qh_start)>(&qh_start),
				"start reduced momentum transfer q_h"));
			args.add(boost::make_shared<opts::option_description>(
				"qk_start", opts::value<decltype(qk_start)>(&qk_start),
				"start reduced momentum transfer q_k"));
			args.add(boost::make_shared<opts::option_description>(
				"ql_start", opts::value<decltype(ql_start)>(&ql_start),
				"start reduced momentum transfer q_l"));

			args.add(boost::make_shared<opts::option_description>(
				"qh_end", opts::value<decltype(qh_end)>(&qh_end),
				"end reduced momentum transfer q_h"));
			args.add(boost::make_shared<opts::option_description>(
				"qk_end", opts::value<decltype(qk_end)>(&qk_end),
				"end reduced momentum transfer q_k"));
			args.add(boost::make_shared<opts::option_description>(
				"ql_end", opts::value<decltype(ql_end)>(&ql_end),
				"end reduced momentum transfer q_l"));

			args.add(boost::make_shared<opts::option_description>(
				"Rx", opts::value<decltype(Rx)>(&Rx),
				"q rotation axis x component"));
			args.add(boost::make_shared<opts::option_description>(
				"Ry", opts::value<decltype(Ry)>(&Ry),
				"q rotation axis y component"));
			args.add(boost::make_shared<opts::option_description>(
				"Rz", opts::value<decltype(Rz)>(&Rz),
				"q rotation axis z component"));
			args.add(boost::make_shared<opts::option_description>(
				"Ralpha", opts::value<decltype(Ralpha)>(&Ralpha),
				"q rotation angle"));

			args.add(boost::make_shared<opts::option_description>(
				"num_threads", opts::value<decltype(num_threads)>(&num_threads),
				"number of threads for calculation"));


			clparser.options(args);
			opts::basic_parsed_options<char> parsedopts = clparser.run();

			opts::variables_map opts_map;
			opts::store(parsedopts, opts_map);
			opts::notify(opts_map);

			if(show_help)
			{
				std::cout << args << std::endl;
				std::cout << "example usage:\n\t"
					<< argv[0] << " --dyntype=h --use_para_perp_calc=0 --outfile=dyn.dat --Gx=1 --Gy=1 --Gz=0 --Bx=1 --By=1 --Bz=0 --Px=1 --Py=-1 --Pz=0 --T=28.5 --B=0.15 --num_points=256 --qh_start=-0.2 --qk_start=-0.2 --ql_start=0 --qh_end=0.2 --qk_end=0.2 --ql_end=0 --Rx=0 --Ry=0 --Rz=1 --Ralpha=0\n"
					<< std::endl;
				return 0;
			}
		}
		catch(const std::exception& ex)
		{
			std::cerr << ex.what() << std::endl;
			return -1;
		}
	}


	if(use_para_perp_calc)
	{
		// fixed number of threads in the old function
		std::cout << "Calculating high-symmetry path dispersion in 4 threads..." << std::endl;

		// calculate the dispersion using simple parallel or perpendicular momentum transfers
		calc_disp_para_perp(dyntype, do_proj,
			Gx,Gy,Gz, Bx,By,Bz, Px,Py,Pz,
			qperpx,qperpy,qperpz,
			qperp, qperp2,
			outfile, !alongqpara,
			T, B, qrange, qdelta,
			explicit_calc,
			gs_file);
	}
	else
	{
		std::cout << "Calculating arbitrary path dispersion in " << num_threads << " threads..." << std::endl;

		// calculate the dispersion using arbitrary momentum transfer paths
		calc_disp_path(dyntype, do_proj,
			Gx,Gy,Gz, Bx,By,Bz, Px,Py,Pz,
			qh_start, qk_start, ql_start,
			qh_end, qk_end, ql_end,
			Rx,Ry,Rz, tl2::d2r(Ralpha),
			num_points, outfile,
			T, B, explicit_calc,
			gs_file, num_threads);
	}

	return 0;
}
