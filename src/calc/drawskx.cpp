/**
 * @author Tobias Weber <tweber@ill.fr>
 * @date mid-2016 (from my phd thesis), oct-2018 (adapted for the paper)
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <fstream>
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
		dAngle = tl2::d2r(dAngle);
		t_mat matRot = tl2::rotation_matrix_3d_z(dAngle);
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


void calcskx(const char* pcFile)
{
	const int iCntX = 24, iCntY = 24;
	const t_real dXScale = 2., dYScale = 2.;
	const t_real dDirScale = 0.3;

	std::ofstream ofstr(pcFile);
	ofstr.precision(8);

	for(int iX=0; iX<iCntX; ++iX)
	{
		t_real dX = -tl2::pi<t_real> + t_real(iX)/t_real(iCntX-1) * 2.*tl2::pi<t_real>;

		for(int iY=0; iY<iCntY; ++iY)
		{
			t_real dY = -tl2::pi<t_real> + t_real(iY)/t_real(iCntY-1) * 2.*tl2::pi<t_real>;

			t_vec vecCoord(3);
			vecCoord[0] = dXScale*dX;
			vecCoord[1] = dYScale*dY;
			vecCoord[2] = 0.;

			t_vec vec = helix_vec(vecCoord)*dDirScale;

			ofstr << vecCoord[0] << "\t" << vecCoord[1] << "\t" << vecCoord[2] << "\t\t";
			ofstr << vec[0] << "\t" << vec[1] << "\t" << vec[2] << std::endl;
		}
	}
}


int main(int argc, char **argv)
{
	calcskx("drawskx.dat");
	return 0;
}
