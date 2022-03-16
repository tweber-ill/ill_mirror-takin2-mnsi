/**
 * Pre-calculates the dispersion and stores it in a database
 * @author Tobias Weber <tweber@ill.fr>
 * @date oct-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "core/skx.h"
#include "tlibs2/libs/str.h"
#include "tlibs2/libs/file.h"
#include <fstream>

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_vec = ublas::vector<t_real>;

#include "core/skx_default_gs.cxx"


// ranges
const t_real qrange = 0.03;
const t_real qstep = 0.0006;
const t_real Erange = 0.5;


void calc_disp(
	t_real qh,
	t_real Bx, t_real By, t_real Bz,
	t_real Px, t_real Py, t_real Pz,
	const std::string& filename)
{
	Skx<t_real, t_cplx, DEF_SKX_ORDER> skx;

	std::vector<ublas::vector<t_cplx>> fourier_skx;
	fourier_skx.reserve(_skxgs_allcomps.size()/3);

	for(std::size_t comp=0; comp<_skxgs_allcomps.size(); comp+=3)
		fourier_skx.push_back(tl2::make_vec<ublas::vector<t_cplx>>({_skxgs_allcomps[comp], _skxgs_allcomps[comp+1], _skxgs_allcomps[comp+2]}));

	skx.SetFourier(fourier_skx);
	skx.SetProjNeutron(false);
	skx.SetT(-1000.);
	skx.SetB(25.);
	skx.GenFullFourier();
	skx.SetFilterZeroWeight(true);

	//t_vec G = tl2::make_vec<t_vec>({ 1, 1, 0 });
	t_vec G = tl2::make_vec<t_vec>({ 0, 0, 0 });
	t_vec Pdir = tl2::make_vec<t_vec>({ Px, Py, Pz });
	t_vec Bdir = tl2::make_vec<t_vec>({ Bx, By, Bz });

	skx.SetCoords(Bdir[0],Bdir[1],Bdir[2], Pdir[0],Pdir[1],Pdir[2]);
	skx.SetG(G[0], G[1], G[2]);



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
			std::tie(Es, std::ignore, wsSF1, wsSF2, wsNSF) = skx.GetDisp(Q[0], Q[1], Q[2], -Erange, Erange);

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

	std::string filename = "skxdyn_" + tl2::var_to_str(idx);
	calc_disp(qh, 0,0,1, 1,-1,0, filename);

	return 0;
}
