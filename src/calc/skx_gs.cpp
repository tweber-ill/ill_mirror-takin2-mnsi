/**
 * Calculate the skyrmion ground state fourier components and free energy
 * @author tweber@ill.fr
 * @date aug-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/skx.h"
#include "core/load_gs.h"

#include <iostream>

#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>
namespace opts = boost::program_options;


#define ORDER DEF_SKX_ORDER

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_vec_cplx = ublas::vector<t_cplx>;
constexpr const t_cplx imag(0, 1);


// initial values
// x and y values are the imaginary components
// z values are the real components
const std::vector<t_real> _allcomps =
{{
	// order 1
	/*
	 0., 0., 10., // (0, 0)
	-5., 5., -5., // (1, 0)
	*/

	// order 4
	/*
	0, -0, 10.661399, // (0, 0)
	-5.8528235, 5.8528235, -5.4024996, // (1, 0)
	0.11529868, -0.11529868, 0.086083531, // (2, 0)
	-0.60556343, 0.60556343, -0.62422482, // (2, 1)
	0.097186182, -0.097186182, 0.080249247, // (3, 0)
	0.20886784, -0.20886784, 0.18114454, // (3, 1)
	0.20886783, -0.20886783, 0.18114454, // (3, 2)
	-0.0095813986, 0.0095813986, -0.010114311, // (4, 0)
	-0.0066107407, 0.0066107407, -0.0024093452, // (4, 1)
	0.0090043808, -0.0090043808, 0.010406047, // (4, 2)
	-0.0066107419, 0.0066107419, -0.0024093447, // (4, 3)
	*/

	// order 6
	/*0, -0, 10.663423,
	-5.8470692, 5.8470692, -5.4072537,
	0.11368408, -0.11368408, 0.087342352,
	-0.60662885, 0.60662885, -0.62704011,
	0.093279895, -0.093279895, 0.091380271,
	0.2033343, -0.2033343, 0.1904102,
	0.20333432, -0.20333432, 0.19041021,
	-0.012533017, 0.012533017, -0.011899267,
	-0.0034953845, 0.0034953845, -0.0023337824,
	0.0080275833, -0.0080275833, 0.0097365074,
	-0.0034953539, 0.0034953539, -0.0023337717,
	-0.00095463806, 0.00095463806, -0.0008123956,
	-0.0051899266, 0.0051899266, -0.0047621906,
	-0.0085790015, 0.0085790015, -0.0075715448,
	-0.0085789957, 0.0085789957, -0.0075715372,
	-0.0051899239, 0.0051899239, -0.004762178,
	0.00033241913, -0.00033241913, 0.00035189196,
	0.00065667985, -0.00065667985, 0.00054426709,
	0.00051040889, -0.00051040889, 0.00034827613,
	0.00035820966, -0.00035820966, 0.00019081221,
	0.0005104236, -0.0005104236, 0.00034828652,
	0.00065668922, -0.00065668922, 0.00054427553,*/

	// order 8
	0, -0, 10.663949,
	-5.8467344, 5.8467344, -5.4073865,
	0.11344794, -0.11344794, 0.087241669,
	-0.60675717, 0.60675717, -0.62723458,
	0.093197032, -0.093197032, 0.091528326,
	0.20347625, -0.20347625, 0.19058935,
	0.20347641, -0.20347641, 0.19058924,
	-0.012632488, 0.012632488, -0.01192958,
	-0.0034826404, 0.0034826404, -0.0023551242,
	0.0082454421, -0.0082454421, 0.0097098496,
	-0.003482748, 0.003482748, -0.002355334,
	-0.00083973334, 0.00083973334, -0.00092484802,
	-0.0050330395, 0.0050330395, -0.0049741336,
	-0.0082261063, 0.0082261063, -0.0079687713,
	-0.0082261422, 0.0082261422, -0.0079688385,
	-0.0050330532, 0.0050330532, -0.0049742042,
	0.00041503849, -0.00041503849, 0.00037718803,
	0.00067345635, -0.00067345635, 0.00058717654,
	0.00049263682, -0.00049263682, 0.00037838644,
	0.00031256166, -0.00031256166, 0.00018562695,
	0.00049265724, -0.00049265724, 0.00037838271,
	0.00067347042, -0.00067347042, 0.00058718046,
	-1.2020331e-05, 1.2020331e-05, -3.1236883e-06,
	8.2439728e-05, -8.2439728e-05, 6.4446271e-05,
	0.00023797057, -0.00023797057, 0.00021101146,
	0.00036771722, -0.00036771722, 0.00033117911,
	0.0003677213, -0.0003677213, 0.0003312039,
	0.00023798387, -0.00023798387, 0.00021104407,
	8.2468067e-05, -8.2468067e-05, 6.4446961e-05,
	1.1573236e-05, -1.1573236e-05, 4.5699916e-06,
	-1.9016436e-06, 1.9016436e-06, -5.7036309e-06,
	-1.1445882e-05, 1.1445882e-05, -1.5860363e-05,
	-6.1572296e-06, 6.1572296e-06, -7.3770892e-06,
	2.0303949e-06, -2.0303949e-06, 2.9350973e-07,
	-6.1253728e-06, 6.1253728e-06, -7.3491997e-06,
	-1.1442694e-05, 1.1442694e-05, -1.5854355e-05,
	-1.9149148e-06, 1.9149148e-06, -5.7176191e-06,
}};


int main(int argc, char **argv)
{
	std::cout << "Skyrmion ground state calculator, use --help for options.\n" << std::endl;

	// argument parser
	opts::basic_command_line_parser<char> clparser(argc, argv);
	opts::options_description args("Program arguments");

	bool show_help = false;
	args.add(boost::make_shared<opts::option_description>(
		"help", opts::bool_switch(&show_help), "show program usage"));

	/*t_real T_exp = 28.5;
	args.add(boost::make_shared<opts::option_description>(
		"T_exp", opts::value<decltype(T_exp)>(&T_exp),
		("experimental temperature, default: " + std::to_string(T_exp)).c_str()));*/

	t_real T_theo = -1000;
	args.add(boost::make_shared<opts::option_description>(
		"T_theo", opts::value<decltype(T_theo)>(&T_theo),
		("theoretical temperature, default: " + std::to_string(T_theo)).c_str()));

	t_real B_theo = 25.052945;
	args.add(boost::make_shared<opts::option_description>(
		"B_theo", opts::value<decltype(B_theo)>(&B_theo),
		("theoretical field, default: " + std::to_string(B_theo)).c_str()));

	bool B_half_BC2 = false;
	args.add(boost::make_shared<opts::option_description>(
		"B_half_BC2", opts::bool_switch(&B_half_BC2), "set B to BC2/2"));

	/*bool B_from_exp = false;
	args.add(boost::make_shared<opts::option_description>(
		"B_from_exp", opts::bool_switch(&B_from_exp), "set B from experimental value"));*/

	std::string gs_file = "";
	args.add(boost::make_shared<opts::option_description>(
		"gsfile", opts::value<decltype(gs_file)>(&gs_file),
		"initial ground state file name"));

	std::string outfile = "skx_gs.bin";
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


	Skx<t_real, t_cplx, ORDER> skx;
	skx.SetDebug(true);

	// load a given initial ground state
	std::vector<t_vec_cplx> fourier;
	if(gs_file != "")
	{
		bool ok = false;
		t_real T_theo_file, B_theo_file;
		std::tie(ok, T_theo_file, B_theo_file, fourier, std::ignore) =
			load_gs<std::decay_t<decltype(fourier)>>(gs_file, 's');
		if(!ok)
		{
			std::cerr << "Error: Could not load skx ground state \""
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

	// set default initial fourier components
	else
	{
		fourier.reserve(_allcomps.size()/3);

		for(std::size_t comp=0; comp<_allcomps.size(); comp+=3)
		{
			fourier.push_back(tl2::make_vec<t_vec_cplx>(
			{
				_allcomps[comp] * imag,
				_allcomps[comp+1] * imag,
				_allcomps[comp+2]
			}));
		}
	}

	skx.SetFourier(fourier);
	skx.SetT(T_theo, false);

	if(B_half_BC2)
		B_theo = skx.GetBC2(false)/2.;

	skx.SetB(B_theo, false);
	//skx.SetB(25., false);

	const auto& peaks60 = skx.GetPeaks60();

	std::cout << "Bc2_exp = " << skx.GetBC2(true) << ", "
		<< "Bc2_theo = " << skx.GetBC2(false) << "."
		<< std::endl;

	std::cout.precision(8);
	std::cout << "Order: " << ORDER << std::endl;
	std::cout << "T_theo = " << T_theo << std::endl;
	std::cout << "B_theo = " << B_theo << std::endl;
	std::cout << "F_start = " << skx.F() << std::endl;
	bool ok = skx.minimise(ORDER, 1,1,0, 0,1,1);
	std::cout << "F_min = " << skx.F() << " (ok: " << std::boolalpha << ok << ")" << std::endl;

	std::cout << "\nFourier components for use in skx_default_gs.cxx:\n";
	const auto& fouriers = skx.GetFourier();
	// iterate peaks
	for(std::size_t fourier_idx=0; fourier_idx<fouriers.size(); ++fourier_idx)
	{
		const auto& fourier = fouriers[fourier_idx];
		std::cout << fourier[0].imag() << ", " << fourier[1].imag() << ", " << fourier[2].real() << ", ";
		if(fourier_idx == 0)
		{
			std::cout << "// (0, 0)";
		}
		else
		{
			const auto& pk = peaks60[fourier_idx-1];
			std::cout << "// (" << pk[0] << ", " << pk[1] << ")";
		}
		std::cout << std::endl;
	}

	// also save the fourier components to a binary file
	if(outfile != "")
	{
		if(!save_gs<decltype(fouriers)>(outfile, 's', T_theo, B_theo, fouriers))
			std::cerr << "\nError: Could not open output file \"" << outfile << "\"." << std::endl;
		else
			std::cout << "\nWrote skyrmion ground state to file \"" << outfile << "\"." << std::endl;
	}

	return 0;
}
