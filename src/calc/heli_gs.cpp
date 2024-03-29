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
using t_vec_cplx = ublas::vector<t_cplx>;
constexpr const auto j = t_cplx(0, 1);


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

	t_real B_exp = 0.171;
	args.add(boost::make_shared<opts::option_description>(
		"B_exp", opts::value<decltype(B_theo)>(&B_theo),
		("experimental field, default: " + std::to_string(B_exp)).c_str()));

	t_real scale0 = 10.;
	args.add(boost::make_shared<opts::option_description>(
		"scale0", opts::value<decltype(scale0)>(&scale0),
		("zero-order scaling factor, default: " + std::to_string(scale0)).c_str()));

	t_real scale = 10.;
	args.add(boost::make_shared<opts::option_description>(
		"scale", opts::value<decltype(scale)>(&scale),
		("scaling factor, default: " + std::to_string(scale)).c_str()));

	std::string gs_file = "";
	args.add(boost::make_shared<opts::option_description>(
		"gsfile", opts::value<decltype(gs_file)>(&gs_file),
		"initial ground state file name"));

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

	// load a given initial ground state
	std::vector<t_vec_cplx> fourier;
	if(gs_file != "")
	{
		bool ok = false;
		t_real T_theo_file, B_theo_file;
		std::tie(ok, T_theo_file, B_theo_file, fourier, std::ignore) =
			load_gs<std::decay_t<decltype(fourier)>>(gs_file, 'h');
		if(!ok)
		{
			std::cerr << "Error: Could not load conical ground state \""
				<< gs_file << "\"." << std::endl;
			return -1;
		}

		// if no explicit temperature or field is given on the command line,
		// use the ones from the file
		if(opts_map.find("T_theo") == opts_map.end())
			T_theo = T_theo_file;
		if(opts_map.find("B_theo") == opts_map.end())
			B_theo = B_theo_file;
	}

	// set a default initial ground state
	else
	{
		fourier = std::vector<t_vec_cplx>{
			tl2::make_vec<t_vec_cplx>({0, 0, scale0}),
			// helical order => Re{M} perp. Im{M}
			tl2::make_vec<t_vec_cplx>({1.+j, 1.-j, 0}) / std::sqrt(2) * scale,
 		};
	}

	heli.SetFourier(fourier);
	heli.SetT(T_theo, false);

	t_real bc2 = heli.GetBC2(false);
	if(B_half_BC2)
		B_theo = bc2/2.;
	else if(B_from_exp)
		B_theo = get_B_theo(T_theo, T_exp, B_exp, !HELI_USE_HOC);
	heli.SetB(B_theo, false);

	// print the table of fourier components
	std::cout.precision(8);
	std::cout << "Order: " << ORDER << std::endl;
	std::cout << "T_theo = " << T_theo << ", T_exp = " << T_exp << std::endl;
	std::cout << "B_theo = " << B_theo << ", BC2_theo = " << bc2
		<< ", fraction = " << B_theo/bc2 << std::endl;
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
