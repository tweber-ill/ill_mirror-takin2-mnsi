/**
 * tests rotations
 * @author tweber@ill.fr
 * @date may-2022
 * @license GPLv2 (see 'LICENSE' file)
 *
 * g++ -std=c++17 -I ../../ext/ -o test_rot test_rot.cpp
 */

#include <iostream>
#include "tlibs2/libs/math17.h"

using t_real = double;
using t_mat = ublas::matrix<t_real>;
using t_vec = ublas::vector<t_real>;
using t_quat = boost::math::quaternion<t_real>;


int main()
{
	t_real eps = 1e-6;

	t_vec B = tl2::make_vec<t_vec>({ 1., 1., 0. });
	t_vec P = tl2::make_vec<t_vec>({ 1., -1., 0. });

	t_quat quatB = tl2::rotation_quat(B, tl2::make_vec<t_vec>({ 0, 0, 1 }));
	t_vec Pin = tl2::quat_vec_prod(quatB, P);
	t_quat quatPin = tl2::rotation_quat(Pin, tl2::make_vec<t_vec>({ 1, 0, 0 }));
	t_quat rotCoord = quatPin * quatB;

	std::cout << "B = " << B << std::endl;
	std::cout << "P = " << P << std::endl;
	std::cout << "rotation = " << rotCoord << std::endl;

	t_vec mB = -B;
	t_vec mP = -P;

	t_vec Brot = tl2::quat_vec_prod(rotCoord, B);
	t_vec Prot = tl2::quat_vec_prod(rotCoord, P);
	t_vec mBrot = tl2::quat_vec_prod(rotCoord, mB);
	t_vec mProt = tl2::quat_vec_prod(rotCoord, mP);

	tl2::set_eps_0(Brot, eps);
	tl2::set_eps_0(Prot, eps);
	tl2::set_eps_0(mBrot, eps);
	tl2::set_eps_0(mProt, eps);

	std::cout << "rotation(B)  = " << Brot << std::endl;
	std::cout << "rotation(P)  = " << Prot << std::endl;
	std::cout << "rotation(-B) = " << mBrot << std::endl;
	std::cout << "rotation(-P) = " << mProt << std::endl;

	return 0;
}
