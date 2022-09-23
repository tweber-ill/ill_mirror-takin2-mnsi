/**
 * Calculate the conical/helimagnetic ground state fourier components and free energy
 * @author tweber@ill.fr
 * @date aug-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/heli.h"
#include "core/load_gs.h"

#include <iostream>

#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>
namespace opts = boost::program_options;


#define ORDER DEF_HELI_ORDER

using t_real = double;
using t_cplx = std::complex<t_real>;


int main(int argc, char** argv)
{
	std::cout << "Helical/conical ground state calculator, use --help for options.\n" << std::endl;

	// argument parser
	opts::basic_command_line_parser<char> clparser(argc, argv);
	opts::options_description args("Program arguments");

	bool show_help = false;
	args.add(boost::make_shared<opts::option_description>(
		"help", opts::bool_switch(&show_help), "show program usage"));

	t_real T_exp = 28.5;
	args.add(boost::make_shared<opts::option_description>(
		"T_exp", opts::value<decltype(T_exp)>(&T_exp),
		("experimental temperature, default: " + std::to_string(T_exp)).c_str()));

	t_real T_theo = -1000;
	args.add(boost::make_shared<opts::option_description>(
		"T_theo", opts::value<decltype(T_theo)>(&T_theo),
		("theoretical temperature, default: " + std::to_string(T_theo)).c_str()));

	t_real B_theo = 22.5;
	args.add(boost::make_shared<opts::option_description>(
		"B_theo", opts::value<decltype(B_theo)>(&B_theo),
		("theoretical field, default: " + std::to_string(B_theo)).c_str()));

	bool B_half_BC2 = false;
	args.add(boost::make_shared<opts::option_description>(
		"B_half_BC2", opts::bool_switch(&B_half_BC2), "set B to BC2/2"));

	bool B_from_exp = false;
	args.add(boost::make_shared<opts::option_description>(
		"B_from_exp", opts::bool_switch(&B_from_exp), "set B from experimental value"));

	t_real scale0 = 10.;
	args.add(boost::make_shared<opts::option_description>(
		"scale0", opts::value<decltype(scale0)>(&scale0),
		("zero-order scaling factor, default: " + std::to_string(scale0)).c_str()));

	t_real scale = 10.;
	args.add(boost::make_shared<opts::option_description>(
		"scale", opts::value<decltype(scale)>(&scale),
		("scaling factor, default: " + std::to_string(scale)).c_str()));

	std::string outfile = "heli_gs.bin";
	args.add(boost::make_shared<opts::option_description>(
		"outfile", opts::value<decltype(outfile)>(&outfile),
		("output file, default: " + outfile).c_str()));


	clparser.options(args);
	opts::basic_parsed_options<char> parsedopts = clparser.run();

	opts::variables_map opts_map;
	opts::store(parsedopts, opts_map);
	opts::notify(opts_map);

	if(show_help)
	{
		std::cout << args << std::endl;
		return 0;
	}


	Heli<t_real, t_cplx, DEF_HELI_ORDER> heli;

	const auto j = t_cplx(0, 1);
	std::vector<ublas::vector<t_cplx>> fourier{
		tl2::make_vec<ublas::vector<t_cplx>>({0, 0, scale0}),
		// helical order => Re{M} perp. Im{M}
		tl2::make_vec<ublas::vector<t_cplx>>({1.+j, 1.-j, 0}) / std::sqrt(2) * scale,
 	};

	heli.SetFourier(fourier);
	heli.SetT(T_theo, false);

	if(B_half_BC2)
		B_theo = heli.GetBC2(false)/2.;
	else if(B_from_exp)
		B_theo = get_B_theo(T_theo, T_exp, 0.171, !HELI_USE_HOC);

	heli.SetB(B_theo, false);

	// print the table of fourier components
	std::cout.precision(8);
	std::cout << "Order: " << ORDER << std::endl;
	std::cout << "T_theo = " << T_theo << ", T_exp = " << T_exp << std::endl;
	std::cout << "B_theo = " << B_theo << std::endl;
	std::cout << "F_start = " << heli.F() << std::endl;
	bool ok = heli.minimise(ORDER, 0,1,0, 0,1,1);
	std::cout << "F_min = " << heli.F() << " (ok: " << std::boolalpha << ok << ")" << std::endl;

	std::cout << "\nFourier components for use in heli_default_gs.cxx:\n";
	const auto& fouriers = heli.GetFourier();
	// iterate peaks
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

	// also save the fourier components to a binary file
	if(outfile != "")
	{
		if(!save_gs<decltype(fouriers)>(outfile, 'h', T_theo, B_theo, fouriers))
			std::cerr << "\nError: Could not open output file \"" << outfile << "\"." << std::endl;
		else
			std::cout << "\nWrote conical ground state to file \"" << outfile << "\"." << std::endl;
	}

	return 0;
}
