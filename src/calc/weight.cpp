/**
 * Calculates the weights/energies at a specific q
 * @author Tobias Weber <tweber@ill.fr>
 * @date sep-19
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/heli.h"
#include "core/skx.h"
#include "core/fp.h"
#include "core/load_gs.h"

#include <fstream>
#include <memory>

#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>
namespace opts = boost::program_options;

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_vec = ublas::vector<t_real>;
using t_vec_cplx = ublas::vector<t_cplx>;

#include "core/skx_default_gs.cxx"
#include "core/heli_default_gs.cxx"


static bool calc_weight(char dyntype,
	t_real Gx, t_real Gy, t_real Gz,
	t_real Bx, t_real By, t_real Bz,
	t_real Px, t_real Py, t_real Pz,
	t_real Qx, t_real Qy, t_real Qz,
	bool do_proj = true, t_real T = -1., t_real B = -1.,
	const std::string& gs_file = "", bool explicit_calc = true,
	std::ostream& ostr = std::cout)
{
	std::shared_ptr<MagDynamics<t_real, t_cplx>> dyn;

	if(dyntype == 's')
	{
		// get default ground state
		auto [skxgs_T, skxgs_B, skxgs] = _get_skx_gs<t_vec_cplx>();

		// optionally load a given ground state file
		if(gs_file != "")
		{
			bool ok = false;
			std::tie(ok, skxgs_T, skxgs_B, skxgs, std::ignore) =
				load_gs<std::decay_t<decltype(skxgs)>>(gs_file, 's');
			if(!ok)
			{
				std::cerr << "Error: Could not load skyrmion ground state \""
					<< gs_file << "\"." << std::endl;
				return false;
			}
		}

		auto skx = std::make_shared<Skx<t_real, t_cplx, DEF_SKX_ORDER>>();
		skx->SetFourier(skxgs);
		skx->SetT(skxgs_T, false);
		skx->SetB(/*skx->GetBC2(false)/2.*/ skxgs_B, false);

		dyn = skx;
	}
	else if(dyntype == 'h')
	{
		// get default ground state
		auto [heligs_T, heligs_B, heligs] = _get_heli_gs<t_vec_cplx>();

		// optionally load a given ground state file
		if(gs_file != "")
		{
			bool ok = false;
			std::tie(ok, heligs_T, heligs_B, heligs, std::ignore) =
				load_gs<std::decay_t<decltype(heligs)>>(gs_file, 'h');
			if(!ok)
			{
				std::cerr << "Error: Could not load conical ground state \""
					<< gs_file << "\"." << std::endl;
				return false;
			}
		}

		auto heli = std::make_shared<Heli<t_real, t_cplx, DEF_HELI_ORDER>>();
		heli->SetExplicitCalc(explicit_calc);
		heli->SetFourier(heligs);
		heli->SetT(heligs_T, false);
		heli->SetB(/*heli->GetBC2(false)/2.*/ heligs_B, false);

		dyn = heli;
	}
	else if(dyntype == 'f')
	{
		dyn = std::make_shared<FP<t_real, t_cplx>>();
	}
	else
	{
		std::cerr << "Unknown dynamics type selected." << std::endl;
		return false;
	}


	dyn->SetCoords(Bx,By,Bz, Px,Py,Pz);
	dyn->SetT(T, true);
	dyn->SetB(B, true);
	dyn->SetFilterZeroWeight(1);
	dyn->SetProjNeutron(do_proj);
	dyn->SetG(Gx, Gy, Gz);


	t_real Erange = -1.;	// negative: disable range
	auto [Es, wsUnpol, wsSF1, wsSF2, wsNSF] = dyn->GetDisp(Qx, Qy, Qz, -Erange, Erange);

	ostr
		<< "# B = " << B << "\n"
		<< "# T = " << T << "\n"
		<< "# Magnetic phase: " << dyntype << "\n"
		<< std::setw(15) << std::left << "# No." << " "
		<< std::setw(15) << std::left << "E (meV)" << " "
		<< std::setw(15) << std::left << "w_total" << " "
		<< std::setw(15) << std::left << "w_SF1" << " "
		<< std::setw(15) << std::left << "w_SF2" << " "
		<< std::setw(15) << std::left << "w_NSF" << "\n";

	for(std::size_t i=0; i<Es.size(); ++i)
	{
		ostr
			<< std::setw(15) << std::left << (i+1) << " "
			<< std::setw(15) << std::left /*<< std::scientific*/ << Es[i] << " "
			<< std::setw(15) << std::left /*<< std::scientific*/ << wsUnpol[i] << " "
			<< std::setw(15) << std::left /*<< std::scientific*/ << wsSF1[i] << " "
			<< std::setw(15) << std::left /*<< std::scientific*/ << wsSF2[i] << " "
			<< std::setw(15) << std::left /*<< std::scientific*/ << wsNSF[i] << "\n";
	}

	return true;
}


int main(int argc, char** argv)
{
	std::cout
		<< "--------------------------------------------------------------------------------\n"
		<< "\tDynamical structure factor calculation tool,\n\t\tT. Weber <tweber@ill.fr>, September 2019.\n"
		<< "--------------------------------------------------------------------------------\n\n";

	t_real Gx = 1., Gy = 1., Gz = 0.;
	t_real Bx = 1., By = 1., Bz = 0.;
	t_real Qx = 1.068, Qy = 0.932, Qz = 0.;
	t_real Px = 1., Py = -1., Pz = 0.;
	t_real B = 0.17, T = 20.;

	char dyntype = 'h';
	std::string gs_file = "";
	std::string outfile = "";

	bool explicit_calc = true;
	bool do_proj = true;

	if(argc <= 1)
	{
		std::cout << "No arguments given, running interactively. "
			<< "Use \"" << argv[0] << " --help\" to show all options.\n"
			<< std::endl;

		while(true)
		{
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
			std::cin >> do_proj;

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
			std::cout << "# Q_proj = " << do_proj << "\n";
		}
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
				"dyntype", opts::value<decltype(dyntype)>(&dyntype),
				"dispersion type [s/h/f]"));

			args.add(boost::make_shared<opts::option_description>(
				"explicit_calc", opts::value<decltype(explicit_calc)>(&explicit_calc),
				"use explicit calculation"));

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
				"Qx", opts::value<decltype(Qx)>(&Qx),
				"momentum transfer vector x component"));
			args.add(boost::make_shared<opts::option_description>(
				"Qy", opts::value<decltype(Qy)>(&Qy),
				"momentum transfer vector y component"));
			args.add(boost::make_shared<opts::option_description>(
				"Qz", opts::value<decltype(Qz)>(&Qz),
				"momentum transfer vector z component"));

			args.add(boost::make_shared<opts::option_description>(
				"T", opts::value<decltype(T)>(&T),
				"temperature"));
			args.add(boost::make_shared<opts::option_description>(
				"B", opts::value<decltype(B)>(&B),
				"magnetic field magnitude"));


			clparser.options(args);
			opts::basic_parsed_options<char> parsedopts = clparser.run();

			opts::variables_map opts_map;
			opts::store(parsedopts, opts_map);
			opts::notify(opts_map);

			if(show_help)
			{
				std::cout << args << std::endl;
				std::cout << "example usage:\n\t"
					<< argv[0] << "  --dyntype=h --outfile=weight.dat --Gx=1 --Gy=1 --Gz=0 --Bx=1 --By=1 --Bz=0 --Px=1 --Py=-1 --Pz=0 --Qx=1.05 --Qy=1.05 --Qz=0 --T=28.5 --B=0.15\n"
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

	// optional output file
	std::ostream *ostr = &std::cout;
	std::ofstream ofstr;
	if(outfile != "")
	{
		ofstr = std::ofstream(outfile);
		ostr = &ofstr;

		std::cout << "Writing results to \"" << outfile << "\"."
			<< std::endl;
	}
	ostr->precision(8);

	if(!calc_weight(dyntype, Gx,Gy,Gz, Bx,By,Bz, Px,Py,Pz, Qx,Qy,Qz,
		do_proj, T, B, gs_file, explicit_calc, *ostr))
	{
		std::cerr << "Calculation failed!" << std::endl;
		return -1;
	}

	std::cout << std::endl;
	return 0;
}
