/**
 * loads and saves ground states
 * @author Tobias Weber <tweber@ill.fr>
 * @date 23/sep/2022
 * @license GPLv2 (see 'LICENSE' file)
 */

#ifndef __LOAD_GS_H__
#define __LOAD_GS_H__

#include <iostream>
#include <fstream>
#include <tuple>

#include "tlibs2/libs/math17.h"


/**
 * load ground state configuration
 */
template<class t_fouriers, class t_real = std::decay_t<decltype(t_fouriers{}[0][0].real())>>
std::tuple<bool, t_real /*T*/, t_real /*B*/, t_fouriers, char /*gs_type*/>
load_gs(const std::string& infile, char expected_gs_type = '*')
{
	std::ifstream ifstr(infile, std::ios_base::binary);
	if(!ifstr)
	{
		std::cerr << "Error: Could not open input ground state file \""
			<< infile << "\"." << std::endl;
		return std::make_tuple(false, 0., 0., t_fouriers{}, '-');
	}

	// check magic values
	char magic[4];
	ifstr.read(magic, sizeof(magic));
	if(magic[0]!='g' || magic[1]!='s' || magic[2]!='_')
	{
		std::cerr << "Error: File \"" << infile << "\""
			<< " does not have a valid ground state." << std::endl;
		return std::make_tuple(false, 0., 0., t_fouriers{}, '-');
	}

	// check ground state type
	if(expected_gs_type != '*' && expected_gs_type != magic[3])
	{
		std::cerr << "Error: File \"" << infile << "\""
			<< " has a mismatching ground state type." << std::endl;
		return std::make_tuple(false, 0., 0., t_fouriers{}, magic[3]);
	}

	// read temperature and field
	t_real T_theo = 0.;
	t_real B_theo = 0.;
	ifstr.read((char*)&T_theo, sizeof(T_theo));
	ifstr.read((char*)&B_theo, sizeof(B_theo));

	// read fourier components
	t_fouriers fouriers;

	using t_vec_cplx = std::decay_t<decltype(fouriers[0])>;
	using t_cplx = std::decay_t<decltype(t_vec_cplx{}[0])>;

	while(true)
	{
		t_real components[6]{};
		if(!ifstr.read((char*)components, sizeof(components)))
			break;

		fouriers.emplace_back(tl2::make_vec<t_vec_cplx>(
		{
			t_cplx{components[0], components[1]},
			t_cplx{components[2], components[3]},
			t_cplx{components[4], components[5]}
		}));
	}

	return std::make_tuple(true, T_theo, B_theo, fouriers, magic[3]);
}


/**
 * save ground state configuration
 */
template<class t_fouriers, class t_real = std::decay_t<decltype(t_fouriers{}[0][0].real())>>
bool save_gs(const std::string& outfile, char gs_type,
	t_real T_theo, t_real B_theo, const t_fouriers& fouriers)
{
	std::ofstream ofstr(outfile, std::ios_base::binary);
	if(!ofstr)
	{
		std::cerr << "Error: Could not open output ground state file \""
			<< outfile << "\"." << std::endl;
		return false;
	}

	// write magic values
	const char magic[] = {'g', 's', '_', gs_type};
	ofstr.write(magic, sizeof(magic));

	// write temperature and field
	ofstr.write((char*)&T_theo, sizeof(T_theo));
	ofstr.write((char*)&B_theo, sizeof(B_theo));

	// iterate peaks
	for(std::size_t pk_idx=0; pk_idx<fouriers.size(); ++pk_idx)
	{
		// write fourier components
		t_real x_re = fouriers[pk_idx][0].real();
		t_real x_im = fouriers[pk_idx][0].imag();
		t_real y_re = fouriers[pk_idx][1].real();
		t_real y_im = fouriers[pk_idx][1].imag();
		t_real z_re = fouriers[pk_idx][2].real();
		t_real z_im = fouriers[pk_idx][2].imag();

		ofstr.write((char*)&x_re, sizeof(x_re));
		ofstr.write((char*)&x_im, sizeof(x_im));
		ofstr.write((char*)&y_re, sizeof(y_re));
		ofstr.write((char*)&y_im, sizeof(y_im));
		ofstr.write((char*)&z_re, sizeof(z_re));
		ofstr.write((char*)&z_im, sizeof(z_im));
	}

	return true;
}


#endif
