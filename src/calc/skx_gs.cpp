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
	0, -0, 11.115651, // (0, 0)
	-5.746439, 5.746439, -5.3380396, // (1, 0)
	0.023270633, -0.023270633, 0.0017156114, // (2, 0)
	-0.71767595, 0.71767595, -0.72431416, // (2, 1)
	0.10929854, -0.10929854, 0.10703773, // (3, 0)
	0.20599633, -0.20599633, 0.19329563, // (3, 1)
	0.20599633, -0.20599633, 0.19329556, // (3, 2)
	-0.011716128, 0.011716128, -0.011024157, // (4, 0)
	0.0040674794, -0.0040674794, 0.0049866219, // (4, 1)
	0.018732234, -0.018732234, 0.019801169, // (4, 2)
	0.0040674151, -0.0040674151, 0.0049865138, // (4, 3)
	-0.0019807787, 0.0019807787, -0.0019822352, // (5, 0)
	-0.0069159196, 0.0069159196, -0.0068085202, // (5, 1)
	-0.0099036115, 0.0099036115, -0.0096291809, // (5, 2)
	-0.0099036453, 0.0099036453, -0.0096292126, // (5, 3)
	-0.0069159607, 0.0069159607, -0.006808545, // (5, 4)
	0.00060021326, -0.00060021326, 0.00060406467, // (6, 0)
	0.00061130456, -0.00061130456, 0.00059018794, // (6, 1)
	-1.3735676e-05, 1.3735676e-05, -6.6229505e-05, // (6, 2)
	-0.00042105525, 0.00042105525, -0.00047811525, // (6, 3)
	-1.3726923e-05, 1.3726923e-05, -6.6221851e-05, // (6, 4)
	0.00061131503, -0.00061131503, 0.00059019741, // (6, 5)
	-1.3980691e-05, 1.3980691e-05, -1.7719541e-05, // (7, 0)
	0.00017636166, -0.00017636166, 0.00016235658, // (7, 1)
	0.00044079235, -0.00044079235, 0.00039095892, // (7, 2)
	0.00060436339, -0.00060436339, 0.0005279182, // (7, 3)
	0.00060435248, -0.00060435248, 0.00052792323, // (7, 4)
	0.00044079392, -0.00044079392, 0.00039096036, // (7, 5)
	0.00017635893, -0.00017635893, 0.00016235066, // (7, 6)
	-2.1648297e-05, 2.1648297e-05, -2.3130656e-05, // (8, 0)
	-5.292669e-05, 5.292669e-05, -4.4138977e-05, // (8, 1)
	-5.0631786e-05, 5.0631786e-05, -3.760949e-05, // (8, 2)
	-2.1828066e-05, 2.1828066e-05, -1.0035428e-05, // (8, 3)
	-3.2241895e-06, 3.2241895e-06, 5.8179637e-06, // (8, 4)
	-2.1822282e-05, 2.1822282e-05, -1.0031356e-05, // (8, 5)
	-5.0638447e-05, 5.0638447e-05, -3.7612525e-05, // (8, 6)
	-5.2927678e-05, 5.2927678e-05, -4.4142297e-05, // (8, 7)
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

	t_real B_theo = 23.9283; //25.052945;
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
