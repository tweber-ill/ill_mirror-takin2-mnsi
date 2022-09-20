/**
 * converts version 1 grid data to version 2
 * @author Tobias Weber <tweber@ill.fr>
 * @date 05-jan-20
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>
namespace opts = boost::program_options;


using t_float_src = double;  // input grid data type
using t_float_dst = double;  // output grid data type


static bool convert(
	const std::string& filenameIdx, const std::string& filenameDat, // input files
	const std::string& filenameNewDat,                              // output files
	t_float_src eps, t_float_src maxE, int channel,
	t_float_dst hstep, t_float_dst hmin, t_float_dst hmax,          //
	t_float_dst kstep, t_float_dst kmin, t_float_dst kmax,          // input file dimensions
	t_float_dst lstep, t_float_dst lmin, t_float_dst lmax)          //
{
	std::cout << "Converting data file \"" << filenameDat << "\"..." << std::endl;

	std::ifstream ifDat(filenameDat);
	if(!ifDat)
	{
		std::cerr << "Error: Cannot open data file \"" << filenameDat << "\"." << std::endl;
		return false;
	}

	std::ofstream ofDat(filenameNewDat);
	if(!ofDat)
	{
		std::cerr << "Error: Cannot open output file \"" << filenameNewDat << "\"." << std::endl;
		return false;
	}

	std::unordered_map<std::size_t, std::size_t> mapIndices;
	std::size_t removedBranches = 0;


	// write a dummy index file offset at the beginning (to be filled in later)
	std::size_t idx_offs = 0;
	ofDat.write((char*)&idx_offs, sizeof(idx_offs));

	// write data dimensions
	ofDat.write((char*)&hmin, sizeof(hmin));
	ofDat.write((char*)&hmax, sizeof(hmax));
	ofDat.write((char*)&hstep, sizeof(hstep));

	ofDat.write((char*)&kmin, sizeof(kmin));
	ofDat.write((char*)&kmax, sizeof(kmax));
	ofDat.write((char*)&kstep, sizeof(kstep));

	ofDat.write((char*)&lmin, sizeof(lmin));
	ofDat.write((char*)&lmax, sizeof(lmax));
	ofDat.write((char*)&lstep, sizeof(lstep));

	// header
	ofDat << "takin_grid_data_ver2|title:skxdyn|author:tweber@ill.fr|date:25/mar/2020";

	ifDat.seekg(0, std::ios_base::end);
	std::streampos ifDat_size = ifDat.tellg();
	ifDat.seekg(0, std::ios_base::beg);

	std::size_t iLoop = 0;
	int last_progress = -1;

	while(1)
	{
		std::size_t idx = ifDat.tellg();
		std::size_t idxnew = ofDat.tellp();

		// show progress
		int cur_progress = (idx*100) / ifDat_size;
		if(cur_progress != last_progress)
		{
			std::cout << "Running: " << cur_progress << " % -- "
				<< idx/1024/1024 << " MB read, "
				<< idxnew/1024/1024 << " MB written."
				<< "            \r"
				<< std::flush;
			last_progress = cur_progress;
		}

		unsigned int numBranches = 0;
		ifDat.read((char*)&numBranches, sizeof(numBranches));
		if(ifDat.gcount() != sizeof(numBranches) || ifDat.eof())
			break;

		mapIndices.insert(std::make_pair(idx, idxnew));

		// write placeholder
		unsigned int numNewBranches = 0;
		ofDat.write((char*)&numNewBranches, sizeof(numNewBranches));


		for(unsigned int branch=0; branch<numBranches; ++branch)
		{
			t_float_src vals[4] = { 0, 0, 0, 0 };
			ifDat.read((char*)vals, sizeof(vals));

			t_float_dst E = t_float_dst(vals[0]);
			if(std::abs(E) < eps)
				E = 0.;

			t_float_dst w[3] =
			{
				t_float_dst(std::abs(vals[1])),
				t_float_dst(std::abs(vals[2])),
				t_float_dst(std::abs(vals[3])),
			};
			t_float_dst total_w = w[0] + w[1] + w[2];

			// do some operation on the channels
			if(channel < 0)
			{
				if(total_w >= eps && std::abs(E) <= maxE)
				{
					if(channel == -1)       // keep original channels
					{
						t_float_dst newvals[4] = { E, w[0], w[1], w[2] };
						ofDat.write((char*)newvals, sizeof(newvals));
					}
					else if(channel == -2)  // add all channels
					{
						t_float_dst newvals[2] = { E, total_w };
						ofDat.write((char*)newvals, sizeof(newvals));
					}

					++numNewBranches;
				}
				else
				{
					++removedBranches;
				}
			}

			// pick a single channel
			else if(channel >= 0)
			{
				if(w[channel] >= eps && std::abs(E) <= maxE)
				{
					t_float_dst newvals[2] = { E, w[channel] };
					ofDat.write((char*)newvals, sizeof(newvals));

					++numNewBranches;
				}
				else
				{
					++removedBranches;
				}
			}
		}

		// seek back and write real number of branches
		std::size_t lastIdx = ofDat.tellp();
		ofDat.seekp(idxnew, std::ios_base::beg);
		ofDat.write((char*)&numNewBranches, sizeof(numNewBranches));
		ofDat.seekp(lastIdx, std::ios_base::beg);

		++iLoop;
	}

	ifDat.close();

	// update index at beginning
	idx_offs = ofDat.tellp();
	ofDat.seekp(0, std::ios_base::beg);
	ofDat.write((char*)&idx_offs, sizeof(idx_offs));
	ofDat.seekp(idx_offs, std::ios_base::beg);

	std::cout << removedBranches << " branches removed (w<eps || E<maxE)." << std::endl;
	std::cout << "\nConverting index file \"" << filenameIdx << "\"..." << std::endl;

	std::ifstream ifIdx(filenameIdx);
	if(!ifIdx)
	{
		std::cerr << "Error: Cannot open index file \"" << filenameIdx << "\"." << std::endl;
		return false;
	}

	std::ofstream &ofIdx = ofDat;

	while(1)
	{
		std::size_t idx = 0;

		ifIdx.read((char*)&idx, sizeof(idx));
		if(ifIdx.gcount() != sizeof(idx) || ifIdx.eof())
			break;

		auto iter = mapIndices.find(idx);
		if(iter == mapIndices.end())
		{
			std::cerr << "Error: Index " << std::hex << idx << " was not found." << std::endl;
			continue;
		}

		std::size_t newidx = iter->second;
		ofIdx.write((char*)&newidx, sizeof(newidx));
	}

	return true;
}


int main(int argc, char** argv)
{
	std::cout << "Takin grid version 1 to version 2 converter.\n" << std::endl;

	// arguments
	bool show_help = false;

	std::string filenameIdx = "skxdyn.idx";
	std::string filenameDat = "skxdyn.bin";
	std::string filenameNewDat = "skxdyn.sqw";

	t_float_src eps = 1e-8;
	t_float_src maxE = 1.5;
	int channel = -1;

	t_float_dst hstep = 0.0006;
	t_float_dst hmin = -0.03;
	t_float_dst hmax = 0.03;

	t_float_dst kstep = 0.0006;
	t_float_dst kmin = -0.03;
	t_float_dst kmax = 0.03;

	t_float_dst lstep = 0.0006;
	t_float_dst lmin = -0.03;
	t_float_dst lmax = 0.03;

	// argument parser
	opts::basic_command_line_parser<char> clparser(argc, argv);
	opts::options_description args("Program arguments");

	args.add(boost::make_shared<opts::option_description>(
		"help", opts::bool_switch(&show_help), "show program usage"));

	args.add(boost::make_shared<opts::option_description>(
		"infile_idx", opts::value<decltype(filenameIdx)>(&filenameIdx),
		("input index file, default: " + filenameIdx).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"infile_dat", opts::value<decltype(filenameDat)>(&filenameDat),
		("input data file, default: " + filenameDat).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"outfile", opts::value<decltype(filenameNewDat)>(&filenameNewDat),
		("output grid file, default: " + filenameNewDat).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"eps", opts::value<decltype(eps)>(&eps),
		("epsilon value, default: " + std::to_string(eps)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"max_E", opts::value<decltype(maxE)>(&maxE),
		("maximum energy in meV, default: " + std::to_string(maxE)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"channel", opts::value<decltype(channel)>(&channel),
		("polarisation channel (0..2: pick a channel, -1: keep individual channels, -2: add all channels), default: " + std::to_string(channel)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"h_step", opts::value<decltype(hstep)>(&hstep),
		("h step width in rlu, default: " + std::to_string(hstep)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"h_min", opts::value<decltype(hmin)>(&hmin),
		("minimum h in rlu, default: " + std::to_string(hmin)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"h_max", opts::value<decltype(hmax)>(&hmax),
		("maximum h in rlu, default: " + std::to_string(hmax)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"k_step", opts::value<decltype(kstep)>(&kstep),
		("k step width in rlu, default: " + std::to_string(kstep)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"k_min", opts::value<decltype(kmin)>(&kmin),
		("minimum k in rlu, default: " + std::to_string(kmin)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"k_max", opts::value<decltype(kmax)>(&kmax),
		("maximum k in rlu, default: " + std::to_string(kmax)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"l_step", opts::value<decltype(lstep)>(&lstep),
		("l step width in rlu, default: " + std::to_string(lstep)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"l_min", opts::value<decltype(lmin)>(&lmin),
		("minimum l in rlu, default: " + std::to_string(lmin)).c_str()));
	args.add(boost::make_shared<opts::option_description>(
		"l_max", opts::value<decltype(lmax)>(&lmax),
		("maximum l in rlu, default: " + std::to_string(lmax)).c_str()));

	clparser.options(args);
	opts::basic_parsed_options<char> parsedopts = clparser.run();

	opts::variables_map opts_map;
	opts::store(parsedopts, opts_map);
	opts::notify(opts_map);

	if(show_help)
	{
		std::cout << args << std::endl;
		std::cout << "Example usage: " << argv[0]
			<< " --infile_idx=skxdyn_asym.idx --infile_dat=skxdyn_asym.bin --outfile=skxdyn.sqw"
			<< " --h_step=0.001 --k_step=0.001 --l_step=0.001"
			<< " --h_min=-0.096 --k_min=-0.096 --l_min=-0.096"
			<< " --h_max=0.096 --k_max=0.096 --l_max=0.096"
			<< " --channel=0 --eps=2e-2\n" << std::endl;
		return 0;
	}

	if(channel < -2 || channel > 2)
	{
		std::cerr << "Error: Invalid polarisation channel selected." << std::endl;
		return -1;
	}

	// padding
	hmax += hstep;
	kmax += kstep;
	lmax += lstep;

	// start conversion
	bool ok = convert(filenameIdx, filenameDat, filenameNewDat,
		eps, maxE, channel,
		hstep, hmin, hmax,
		kstep, kmin, kmax,
		lstep, lmin, lmax);

	return ok ? 0 : -1;
}
