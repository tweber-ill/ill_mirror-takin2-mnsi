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


template<class t_fouriers, class t_real = decltype(t_fouriers{}[0][0].real())>
bool save_gs(const std::string& outfile, char gs_type,
	t_real T_theo, t_real B_theo, const t_fouriers& fouriers)
{
	std::ofstream ofstr(outfile, std::ios_base::binary);
	if(!ofstr)
	{
		std::cerr << "Error: Could not open output file \"" << outfile << "\"." << std::endl;
		return false;
	}

	// write magic number
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
