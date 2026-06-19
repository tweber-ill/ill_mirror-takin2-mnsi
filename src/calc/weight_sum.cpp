/**
 * calculates the integrated weights/energies for setup 3
 * @author Tobias Weber <tweber@ill.fr>
 * @date jun-20
 * @license GPLv2 (see 'LICENSE' file)
 *
 * TODO: take into account the entire sector, not just the central path
 */

#include "core/heli.h"
#include "core/skx.h"
#include "core/load_gs.h"

#include "tlibs2/libs/phys.h"

#include <fstream>

#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>
#include <boost/histogram.hpp>
namespace opts = boost::program_options;
namespace hist = boost::histogram;

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_vec = ublas::vector<t_real>;
using t_vec_cplx = ublas::vector<t_cplx>;

#include "core/skx_default_gs.cxx"
#include "core/heli_default_gs.cxx"


const t_real g_eps = 1e-5;
const t_real g_weight_eps = 1e-6;

#define E_BINS 200
#define COL_SIZE 15


/**
 * print ground state magnetisation
 */
static void print_groundstate(const std::string& phase, const std::vector<t_vec_cplx>& gs)
{
	std::cout << phase << " ground state:\n";

	for(const auto& fourier : gs)
	{
		std::cout
			<< "{ " << fourier[0].real() << ", " << fourier[0].imag() << " }, "
			<< "{ " << fourier[1].real() << ", " << fourier[1].imag() << " }, "
			<< "{ " << fourier[2].real() << ", " << fourier[2].imag() << " },"
			<< std::endl;
	}

	std::cout << std::endl;
}


void calc_disp(const t_vec& Gvec,
	const t_vec& Bvec, const t_vec& Pvec,
	t_real q, t_real q_oop,
	bool bProj = true, bool filter_zero_weight = true,
	t_real T = 28.5, t_real B = 0.158,
	const std::string& skx_gs_file = "", const std::string& heli_gs_file = "",
	bool explicit_calc = true,
	t_real angle_begin_deg = 0., t_real angle_end_deg = 360., int num_angles = 512,
	std::string outfile = "weightsum")
{
	// vector perpendicular to the pinning, but in the skyrmion plane
	t_vec Pperpvec = tl2::cross_3(Pvec, Bvec);
	Pperpvec /= tl2::veclen(Pperpvec);

	Skx<t_real, t_cplx, DEF_SKX_ORDER> skx;
	Heli<t_real, t_cplx, DEF_HELI_ORDER> heli;


	// get default ground state
	auto [skxgs_T, skxgs_B, skxgs] = _get_skx_gs<t_vec_cplx>();
	auto [heligs_T, heligs_B, heligs] = _get_heli_gs<t_vec_cplx>();

	// optionally load a given skx ground state file
	if(skx_gs_file != "")
	{
		bool ok = false;
		std::tie(ok, skxgs_T, skxgs_B, skxgs, std::ignore) =
			load_gs<std::decay_t<decltype(skxgs)>>(skx_gs_file, 's');

		if(!ok)
		{
			std::cerr << "Error: Could not load skyrmion ground state \""
				<< skx_gs_file << "\"." << std::endl;
			return;
		}
	}

	// optionally load a given helical ground state file
	if(heli_gs_file != "")
	{
		bool ok = false;
		std::tie(ok, heligs_T, heligs_B, heligs, std::ignore) =
			load_gs<std::decay_t<decltype(heligs)>>(heli_gs_file, 'h');

		if(!ok)
		{
			std::cerr << "Error: Could not load helical ground state \""
				<< heli_gs_file << "\"." << std::endl;
			return;
		}
	}

	skx.SetFourier(skxgs);
	heli.SetFourier(heligs);

	heli.SetExplicitCalc(explicit_calc);

	skx.SetProjNeutron(bProj);
	heli.SetProjNeutron(bProj);

	// use the same temperature for both phases
	skx.SetT(skxgs_T, false);
	heli.SetT(skxgs_T, false);
	heli.SetT(T, true);

	// use the same field for both phases
	//t_real bc2 = skx.GetBC2(false);
	skx.SetB(/*bc2/2.*/ skxgs_B, false);
	heli.SetB(/*bc2/2.*/ skxgs_B, false);
	//heli.SetB(0.17, true);
	heli.SetB(B, true);

	skx.SetFilterZeroWeight(filter_zero_weight);
	heli.SetFilterZeroWeight(filter_zero_weight);

	skx.SetWeightEps(g_weight_eps);
	heli.SetWeightEps(g_weight_eps);

	skx.SetCoords(Bvec[0], Bvec[1], Bvec[2], Pvec[0], Pvec[1], Pvec[2]);
	heli.SetCoords(Bvec[0], Bvec[1], Bvec[2]);

	skx.SetG(Gvec[0], Gvec[1], Gvec[2]);
	heli.SetG(Gvec[0], Gvec[1], Gvec[2]);

	t_real Erange = 0.1;

	t_real angle_begin = angle_begin_deg/180.*M_PI;
	t_real angle_end = angle_end_deg/180.*M_PI;
	t_real angle_delta = 2*M_PI/t_real(num_angles);

	auto histWeightsNSF = hist::make_histogram(hist::axis::regular<t_real>(E_BINS, -Erange, Erange, "E"));
	auto histWeightsSF = hist::make_histogram(hist::axis::regular<t_real>(E_BINS, -Erange, Erange, "E"));

	auto histWeightsHeliNSF = hist::make_histogram(hist::axis::regular<t_real>(E_BINS, -Erange, Erange, "E"));
	auto histWeightsHeliSF = hist::make_histogram(hist::axis::regular<t_real>(E_BINS, -Erange, Erange, "E"));

	if(skx.GetFourier().size())
		print_groundstate("Skx", skx.GetFourier());
	if(heli.GetFourier().size())
		print_groundstate("Heli", heli.GetFourier());


	// result files
	std::ofstream ofstr_raw(outfile + "_skx.dat");
	std::ofstream ofstr_raw_heli(outfile + "_heli.dat");

	// write file header
	for(std::ostream* ostr : {&ofstr_raw, &ofstr_raw_heli})
	{
		ostr->precision(8);
		(*ostr)  << std::left << std::setw(COL_SIZE) << "# angle"
			<< " " << std::left << std::setw(COL_SIZE) << "qh"
			<< " " << std::left << std::setw(COL_SIZE) << "qk"
			<< " " << std::left << std::setw(COL_SIZE) << "ql"
			<< " " << std::left << std::setw(COL_SIZE) << "E"
			<< " " << std::left << std::setw(COL_SIZE) << "wSF1"
			<< " " << std::left << std::setw(COL_SIZE) << "wSF2"
			<< " " << std::left << std::setw(COL_SIZE) << "wNSF"
			<< "\n";
	}


	for(t_real angle = angle_begin; angle < angle_end; angle += angle_delta)
	{
		t_vec qvec = q * (Pvec*std::cos(angle) + Pperpvec*std::sin(angle));  // in skx plane
		qvec += q_oop * Bvec;                                                // out of skx plane
		t_vec Qvec = Gvec + qvec;

		std::cout << "# angle: " << angle/M_PI*180.
			<< ", Q = (" << Qvec[0] << ", " << Qvec[1] << ", " << Qvec[2] << ")"
			<< ", |q| = " << tl2::veclen(qvec) << "."
			<< std::endl;

		{
			// skyrmion dispersion
			auto [Es, wsUnpol, wsSF1, wsSF2, wsNSF] = skx.GetDisp(Qvec[0], Qvec[1], Qvec[2], -Erange, Erange);
			for(std::size_t i = 0; i < Es.size(); ++i)
			{
				histWeightsNSF(Es[i], hist::weight(wsNSF[i]*0.5));
				histWeightsSF(Es[i], hist::weight(wsSF1[i]));

				ofstr_raw << std::left << std::setw(COL_SIZE) << angle
					<< " " << std::left << std::setw(COL_SIZE) << qvec[0]
					<< " " << std::left << std::setw(COL_SIZE) << qvec[1]
					<< " " << std::left << std::setw(COL_SIZE) << qvec[2]
					<< " " << std::left << std::setw(COL_SIZE) << Es[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsSF1[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsSF2[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsNSF[i]
					<< std::endl;
			}
		}

		{
			// conical dispersion
			auto [EsH, wsUnpolH, wsSF1H, wsSF2H, wsNSFH] = heli.GetDisp(Qvec[0], Qvec[1], Qvec[2], -Erange, Erange);
			for(std::size_t i=0; i<EsH.size(); ++i)
			{
				histWeightsHeliNSF(EsH[i], hist::weight(wsNSFH[i]*0.5));
				histWeightsHeliSF(EsH[i], hist::weight(wsSF1H[i]));

				ofstr_raw_heli << std::left << std::setw(COL_SIZE) << angle
					<< " " << std::left << std::setw(COL_SIZE) << qvec[0]
					<< " " << std::left << std::setw(COL_SIZE) << qvec[1]
					<< " " << std::left << std::setw(COL_SIZE) << qvec[2]
					<< " " << std::left << std::setw(COL_SIZE) << EsH[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsSF1H[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsSF2H[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsNSFH[i]
					<< std::endl;
			}
		}
	}


	std::ofstream ofstrBinned(outfile + "_bin.dat");
	ofstrBinned.precision(8);

	ofstrBinned << std::left << std::setw(COL_SIZE) << "# E (meV)"
		<< " " << std::left << std::setw(COL_SIZE) << "bose"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_skx_sf"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_skx_nsf"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_skx_sum"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_skx_sum_bose"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_heli_sf"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_heli_nsf"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_heli_sum"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_heli_sum_bose"
		<< "\n";

	auto iterHeliNSF = hist::indexed(histWeightsHeliNSF).begin();
	auto iterHeliSF = hist::indexed(histWeightsHeliSF).begin();
	auto iterNSF = hist::indexed(histWeightsNSF).begin();
	for(const auto& val : hist::indexed(histWeightsSF))
	{
		t_real E = val.bin().lower() + 0.5*(val.bin().upper() - val.bin().lower());
		t_real w = *val / t_real{E_BINS};

		t_real E_nsf = iterNSF->bin().lower() + 0.5*(iterNSF->bin().upper() - iterNSF->bin().lower());
		t_real w_nsf = **iterNSF / t_real{E_BINS};

		t_real E_h = iterHeliSF->bin().lower() + 0.5*(iterHeliSF->bin().upper() - iterHeliSF->bin().lower());
		t_real w_h = **iterHeliSF / t_real{E_BINS};

		t_real E_h_nsf = iterHeliNSF->bin().lower() + 0.5*(iterHeliNSF->bin().upper() - iterHeliNSF->bin().lower());
		t_real w_h_nsf = **iterHeliNSF / t_real{E_BINS};

		if(!tl2::float_equal<t_real>(E, E_h, g_eps) || !tl2::float_equal<t_real>(E, E_h_nsf, g_eps) ||
			!tl2::float_equal<t_real>(E, E_nsf, g_eps))

		{
			std::cerr << "Energy binning mismatch: " << E << " != " << E_h << std::endl;
			break;
		}

		t_real bose = tl2::bose(E, T);

		ofstrBinned << std::left << std::setw(COL_SIZE) << E
			<< " " << std::left << std::setw(COL_SIZE) << bose
			<< " " << std::left << std::setw(COL_SIZE) << w
			<< " " << std::left << std::setw(COL_SIZE) << w_nsf
			<< " " << std::left << std::setw(COL_SIZE) << (w + w_nsf)
			<< " " << std::left << std::setw(COL_SIZE) << (w + w_nsf)*bose
			<< " " << std::left << std::setw(COL_SIZE) << w_h
			<< " " << std::left << std::setw(COL_SIZE) << w_h_nsf
			<< " " << std::left << std::setw(COL_SIZE) << (w_h + w_h_nsf)
			<< " " << std::left << std::setw(COL_SIZE) << (w_h + w_h_nsf)*bose
			<< std::endl;

		++iterHeliSF;
		++iterHeliNSF;
		++iterNSF;
	}
}


int main(int argc, char** argv)
{
	t_real Gx = 0.,  Gy = 0., Gz = 0.;
	t_real Px = 1.,  Py = 1., Pz = 0.;
	t_real Bx = 1., By = -1., Bz = 0.;

	t_real q = 0.0123;    // momentum transfer in skyrmion plane
	t_real q_oop = 0.0;   // momentum transfer out of skyrmion plane

	bool proj = true;     // using the orthogonal 1-|Q><Q| projector
	bool filter_zero_weight = true;
	bool explicit_calc = true;

	std::string skx_gs_file = "";
	std::string heli_gs_file = "";
	std::string outfile = "weightsum";

	t_real T = 28.5;
	t_real B = 0.158;

	// integration arc
	t_real angle_begin = -45;
	t_real angle_end = 270 - 45;
	int num_angles = 512;


	try
	{
		opts::basic_command_line_parser<char> clparser(argc, argv);
		opts::options_description args("program arguments");
		bool show_help = false;

		args.add(boost::make_shared<opts::option_description>(
			"help", opts::bool_switch(&show_help),
			"show program usage"));
		args.add(boost::make_shared<opts::option_description>(
			"explicit_calc", opts::value<decltype(explicit_calc)>(&explicit_calc),
			"use explicit calculation"));
		args.add(boost::make_shared<opts::option_description>(
			"filter_zero_weight", opts::value<decltype(filter_zero_weight)>(&filter_zero_weight),
			"filter out modes with zero spectral weight"));
		args.add(boost::make_shared<opts::option_description>(
			"do_proj", opts::value<decltype(proj)>(&proj),
			"calculate orthogonal projection, S_perp"));
		args.add(boost::make_shared<opts::option_description>(
			"outfile", opts::value<decltype(outfile)>(&outfile),
			"output file name"));
		args.add(boost::make_shared<opts::option_description>(
			"gsfile_skx", opts::value<decltype(skx_gs_file)>(&skx_gs_file),
			"skx ground state file name"));
		args.add(boost::make_shared<opts::option_description>(
			"gsfile_heli", opts::value<decltype(heli_gs_file)>(&heli_gs_file),
			"heli ground state file name"));
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
			"T", opts::value<decltype(T)>(&T),
			"temperature"));
		args.add(boost::make_shared<opts::option_description>(
			"B", opts::value<decltype(B)>(&B),
			"magnetic field magnitude"));
		args.add(boost::make_shared<opts::option_description>(
			"q", opts::value<decltype(q)>(&q),
			"reduced momentum transfer"));
		args.add(boost::make_shared<opts::option_description>(
			"q_oop", opts::value<decltype(q_oop)>(&q_oop),
			"out-of-plane reduced momentum transfer"));
		args.add(boost::make_shared<opts::option_description>(
			"angle_begin", opts::value<decltype(angle_begin)>(&angle_begin),
			"start angle of integration arc"));
		args.add(boost::make_shared<opts::option_description>(
			"angle_end", opts::value<decltype(angle_end)>(&angle_end),
			"end angle of integration arc"));
		args.add(boost::make_shared<opts::option_description>(
			"num_angles", opts::value<decltype(num_angles)>(&num_angles),
			"number of angles for integration arc"));
		//args.add(boost::make_shared<opts::option_description>(
		//	"num_threads", opts::value<decltype(num_threads)>(&num_threads),
		//	"number of threads for calculation"));

		clparser.options(args);
		opts::basic_parsed_options<char> parsedopts = clparser.run();
		opts::variables_map opts_map;
		opts::store(parsedopts, opts_map);
		opts::notify(opts_map);

		if(show_help)
		{
			std::cout << args << std::endl;
			std::cout << "example usage:\n\t"
				<< argv[0] << " --Gx=0 --Gy=0 --Gz=0 --Bx=1 --By=-1 --Bz=0 --Px=1 --Py=1 --Pz=0 --T=28.5 --B=0.15 --num_angles=256 --angle_begin=0 --angle_end=360 --q=0.1\n"
				<< std::endl;
			return 0;
		}
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
		return -1;
	}


	t_vec Gvec = tl2::make_vec<t_vec>({ Gx, Gy, Gz });  // lattice vector
	t_vec Pvec = tl2::make_vec<t_vec>({ Px, Py, Pz });  // pinning vector
	t_vec Bvec = tl2::make_vec<t_vec>({ Bx, By, Bz });  // field direction

	Pvec /= tl2::veclen(Pvec);
	Bvec /= tl2::veclen(Bvec);

	calc_disp(Gvec, Bvec, Pvec,
		q, q_oop, proj, filter_zero_weight,
		T, B, skx_gs_file, heli_gs_file, explicit_calc,
		angle_begin, angle_end, num_angles,
		outfile);
	return 0;
}
