/**
 * @author Tobias Weber <tweber@ill.fr>
 * @date mid-2016 (from my phd thesis), oct-2018 (adapted for the paper)
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <fstream>
#include <vector>
#include <string>

#include "tlibs2/libs/math17.h"


using t_real = double;
using t_vec = ublas::vector<t_real>;
using t_mat = ublas::matrix<t_real>;


t_vec helix_vec(const t_vec& vecCoord)
{
	t_vec vecRet = ublas::zero_vector<t_real>(3);

	// helix angles
	for(t_real dAngle : {0., 120., 240.})
	{
		t_mat matRot = tl2::rotation_matrix_3d_z(tl2::d2r(dAngle));
		t_vec vecCoordRot = ublas::prod(ublas::trans(matRot), vecCoord);

		// helix position
		t_vec vec(3);
		vec[2] = std::cos(vecCoordRot[0]);
		vec[1] = std::sin(vecCoordRot[0]);
		vec[0] = 0;

		vecRet += ublas::prod(matRot, vec);
	}

	vecRet /= ublas::norm_2(vecRet);
	return vecRet;
}


void calcskx(const char* pcFile, int iCnt = 24, t_real dPhaseScale = 2., t_real dVecLen = 0.3)
{
	std::ofstream ofstr(pcFile);
	ofstr.precision(8);

	std::vector<t_real> phi = tl2::linspace<t_real>(
		-tl2::pi<t_real>*dPhaseScale, tl2::pi<t_real>*dPhaseScale, iCnt);

	for(t_real dX : phi)
	{
		for(t_real dY : phi)
		{
			t_vec vecCoord(3);
			vecCoord[0] = dX;
			vecCoord[1] = dY;
			vecCoord[2] = 0.;

			t_vec vec = helix_vec(vecCoord) * dVecLen;

			ofstr << vecCoord[0] << "\t" << vecCoord[1] << "\t" << vecCoord[2] << "\t\t";
			ofstr << vec[0] << "\t" << vec[1] << "\t" << vec[2] << std::endl;
		}
	}
}


int main(int argc, char **argv)
{
	int iCnt = 24;
	t_real dPhaseScale = 2.;
	t_real dVecLen = 0.3;

	if(argc > 1)
		iCnt = std::stoi(argv[1]);
	if(argc > 2)
		dPhaseScale = std::stod(argv[2]);
	if(argc > 3)
		dVecLen = std::stod(argv[3]);

	calcskx("drawskx.dat", iCnt, dPhaseScale, dVecLen);
	return 0;
}
