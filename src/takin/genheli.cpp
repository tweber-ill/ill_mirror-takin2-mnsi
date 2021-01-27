/**
 * Pre-calculates the dispersion and stores it in a database
 * @author Tobias Weber <tweber@ill.fr>
 * @date oct-18, jan-21
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/heli.h"
#include "tlibs2/libs/str.h"
#include "tlibs2/libs/file.h"
#include <fstream>

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_vec = ublas::vector<t_real>;
const auto g_j = t_cplx(0,1);


// ranges
const t_real qrange = 0.1;
const t_real qstep = 0.001;
const t_real Erange = 1.0;


void calc_disp(
	t_real qh,
	t_real Bx, t_real By, t_real Bz,
	const std::string& filename)
{
	Heli<t_real, t_cplx, DEF_HELI_ORDER> heli;

	heli.SetProjNeutron(true);
	heli.SetT(20.);
	heli.SetB(0.3);
	heli.SetFilterZeroWeight(true);

	t_vec G = tl2::make_vec<t_vec>({ 1, 1, 0 });
	t_vec Bdir = tl2::make_vec<t_vec>({ Bx, By, Bz });

	heli.SetCoords(Bdir[0],Bdir[1],Bdir[2]);
	heli.SetG(G[0], G[1], G[2]);



	// outputs
	std::ofstream ofstrLog(filename + ".log");

	if(!ofstrLog)
	{
		std::cerr << "Could not write to output files \"" << filename << ".*\"" << std::endl;
		return;
	}

	ofstrLog << "h = " << qh << std::endl;
	ofstrLog << "k_start = " << -qrange << ", k_end = " << qrange << ", k_step = " << qstep << std::endl;
	ofstrLog << "l_start = " << -qrange << ", l_end = " << qrange << ", l_step = " << qstep << std::endl;


	int k_idx = 0;
	int skip_files = 0;
	for(t_real qk=-qrange; qk<qrange; qk+=qstep)
	{
		std::string thisfile = filename + "_" + tl2::var_to_str(k_idx);
		if(tl2::file_exists<char>((thisfile + ".idx").c_str()) || skip_files > 0)
		{
			std::cout << "Skipping " << thisfile << std::endl;
			++k_idx;
			--skip_files;
			continue;
		}

		std::ofstream ofstrIdx(thisfile + ".idx", std::ios_base::binary);	// format: [h1, k1, l1], idx1,  [h2, k2, l2], idx2, ...
		std::ofstream ofstrBin(thisfile + ".bin", std::ios_base::binary);	// format: num_Es,  E1, w1_SF1, w1_SF2, w1_NSF,  E2, ...


		for(t_real ql=-qrange; ql<qrange; ql+=qstep)
		{
			std::cout << "\rCalculating " << thisfile << " (" << qh << " " << qk << " " << ql << ") ...                ";
			std::cout.flush();

			t_vec Q = G;
			Q[0] += qh; Q[1] += qk; Q[2] += ql;

			std::vector<t_real> Es, wsSF1, wsSF2, wsNSF;
			std::tie(Es, std::ignore, wsSF1, wsSF2, wsNSF) = heli.GetDisp(Q[0], Q[1], Q[2], -Erange, Erange);

			std::size_t offset = (std::size_t)ofstrBin.tellp();
			ofstrIdx.write((char*)&offset, sizeof(offset));

			unsigned int numEs = (unsigned int)Es.size();
			ofstrBin.write((char*)&numEs, sizeof(numEs));

			for(std::size_t idxE=0; idxE<Es.size(); ++idxE)
			{
				ofstrBin.write((char*)&Es[idxE], sizeof(t_real));
				ofstrBin.write((char*)&wsSF1[idxE], sizeof(t_real));
				ofstrBin.write((char*)&wsSF2[idxE], sizeof(t_real));
				ofstrBin.write((char*)&wsNSF[idxE], sizeof(t_real));
			}


			ofstrLog << "q = " << qh << " " << qk <<  " " << ql << ", "
				<< "Q = " << Q[0] << " " << Q[1] <<  " " << Q[2] << ": "
				<< Es.size() << " branches." << std::endl;

			ofstrIdx.flush();
			ofstrBin.flush();
		}

		std::cout << "done" << std::endl;
		++k_idx;
	}
}



int main(int argc, char** argv)
{
	std::cout
		<< "--------------------------------------------------------------------------------\n"
		<< "\tS(q,E) grid calculation tool,\n\t\tby T. Weber <tweber@ill.fr>, October 2018.\n"
		<< "--------------------------------------------------------------------------------\n\n";

	if(argc<=1)
	{
		std::cerr << "Please give an index." << std::endl;
		return -1;
	}

	int idx = tl2::str_to_var<int>(std::string(argv[1]));
	//t_real qh = tl2::str_to_var<t_real>(std::string(argv[2]));
	t_real qh = -qrange + t_real(idx)*qstep;

	std::cout << "qh = " << qh << std::endl;
	std::cout << "idx = " << idx << std::endl;

	std::string filename = "helidyn_" + tl2::var_to_str(idx);
	calc_disp(qh, 1,1,0, filename);

	return 0;
}
