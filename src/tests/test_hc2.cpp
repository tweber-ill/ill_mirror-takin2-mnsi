/**
 * exports hc2 curve
 * @author tweber@ill.fr
 * @date apr-2022
 * @license GPLv2 (see 'LICENSE' file)
 *
 * g++ -std=c++17 -I ../../ext/ -o test_hc2 test_hc2.cpp
 */

#include "../core/constants.h"
#include <iostream>

int main()
{
	double dT = 0.1;

	for(double T = 0.; T < 30.5; T += dT)
	{
		double hc2 = get_bc2<double>(T, false);

		std::cout 
			<< std::left << std::setw(10) << T
			<< std::left << std::setw(10) << hc2
			<< std::endl;
	}

	return 0;
}
