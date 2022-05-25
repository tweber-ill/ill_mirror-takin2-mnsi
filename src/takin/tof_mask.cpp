/**
 * creates a mask for polarisation calulation
 * @author Tobias Weber <tweber@ill.fr>
 * @date may 2022
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <iostream>
#include <string>
#include <cmath>

#include <boost/gil/image.hpp>
#include <boost/gil/extension/io/png.hpp>
namespace gil = boost::gil;

#include "tof.h"
using t_real = double;

//#define USE_BACKGROUND_IMAGE


bool circular_mask(
	const std::string& file,
	t_real rad1, t_real rad2,
	t_real start_angle, t_real end_angle)
{
	gil::gray8_image_t mask(PSD_WIDTH, PSD_HEIGHT);
#ifdef USE_BACKGROUND_IMAGE
	// load a background image to help pose the mask
	gil::read_image("test.png", mask, gil::png_tag());
#endif
	auto mask_view = gil::view(mask);

	t_real x_mid = PSD_WIDTH/2 - 2;
	t_real y_mid = PSD_HEIGHT/2 + 4;

	for(unsigned y=0; y<PSD_HEIGHT; ++y)
	{
		t_real y_pos = y - y_mid;

		auto mask_row = mask_view.row_begin(y);
		for(unsigned x=0; x<PSD_WIDTH; ++x)
		{
			t_real x_pos = x - x_mid;

			t_real rad = std::sqrt(
				x_pos * x_pos +
				y_pos * y_pos);

			t_real angle = std::atan2(y_pos, x_pos);

			bool in_ring = (rad>=rad1 && rad<=rad2);
			bool in_seg = (angle>=start_angle || angle<=end_angle);

#ifdef USE_BACKGROUND_IMAGE
			if(!in_ring || !in_seg)
				continue;
#endif
			mask_row[x] = in_ring && in_seg ? 0xff : 0x00;;
		}
	}

	gil::write_view(file, mask_view, gil::png_tag());
	return true;
}


int main(int argc, char** argv)
{
	bool ok = circular_mask("mask.png",
		34, 46,
		M_PI*0.25, -M_PI*0.25);

	if(ok)
	{
		std::cout << "Created mask file \"mask.png\"." << std::endl;
	}
	else
	{
		std::cerr << "Error writing mask file." << std::endl;
		return -1;
	}

	return 0;
}

