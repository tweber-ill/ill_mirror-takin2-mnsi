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

using t_float_dst = double;


int main()
{
	const double eps = 1e-8;
	const double maxE = 1.5;

	// --------------------------------------------------------------------
	// modify these for each data file to process
	std::string filenameIdx = "/home/tw/tmp/skx/skxdyn.idx";
	std::string filenameDat = "/home/tw/tmp/skx/skxdyn.bin";
	std::string filenameNewDat = "/home/tw/tmp/skx/skxdyn.sqw";

	t_float_dst hstep = 0.0006;
	t_float_dst hmin = -0.03;
	t_float_dst hmax = 0.03 + hstep;

	t_float_dst kstep = 0.0006;
	t_float_dst kmin = -0.03;
	t_float_dst kmax = 0.03 + kstep;

	t_float_dst lstep = 0.0006;
	t_float_dst lmin = -0.03;
	t_float_dst lmax = 0.03 + lstep;
	// --------------------------------------------------------------------


	std::cout << "Converting data file ..." << std::endl;

	std::ifstream ifDat(filenameDat);
	std::ofstream ofDat(filenameNewDat);

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
	ofDat << "takin_grid_data_ver2|title:skxdyn_reseda|author:tweber@ill.fr|date:25/mar/2020";


	std::size_t iLoop = 0;
	while(1)
	{
		std::size_t idx = ifDat.tellg();
		std::size_t idxnew = ofDat.tellp();

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
			double vals[4] = { 0, 0, 0, 0 };
			ifDat.read((char*)vals, sizeof(vals));


			t_float_dst w = t_float_dst(std::abs(vals[1])+std::abs(vals[2])+std::abs(vals[3]));

			if(w >= eps && std::abs(vals[0]) <= maxE)
			{
				t_float_dst E = t_float_dst(vals[0]);
				if(std::abs(E) < eps)
					E = 0.;
				if(std::abs(w) < eps)
					w = 0.;

				t_float_dst newvals[2] = { E, w };
				ofDat.write((char*)newvals, sizeof(newvals));

				++numNewBranches;
			}
			else
			{
				++removedBranches;
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
	std::cout << "\nConverting index file ..." << std::endl;

	std::ifstream ifIdx(filenameIdx);
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

	return 0;
}
