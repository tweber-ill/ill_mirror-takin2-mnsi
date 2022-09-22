/**
 * creates a mask for polarisation calulation
 * @author Tobias Weber <tweber@ill.fr>
 * @date may 2022
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

#include <boost/gil/image.hpp>
#include <boost/gil/extension/io/png.hpp>
namespace gil = boost::gil;

#include "tof.h"
using t_real = double;


bool circular_mask(
	const std::string& file,
	t_real x_mid, t_real y_mid,
	t_real rad1, t_real rad2,
	t_real start_angle, t_real end_angle)
{
	gil::gray8_image_t mask(PSD_WIDTH, PSD_HEIGHT);

	bool has_bckgrd = false;
	fs::path bckgrd_file("test.png");
	if(fs::exists(bckgrd_file))
	{
		// load a background image to help pose the mask
		gil::read_image(bckgrd_file.string(), mask, gil::png_tag());
		has_bckgrd = true;
	}

	auto mask_view = gil::view(mask);

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

			if(has_bckgrd && (!in_ring || !in_seg))
				continue;
			mask_row[x] = in_ring && in_seg ? 0xff : 0x00;;
		}
	}

	gil::write_view(file, mask_view, gil::png_tag());
	return true;
}


int main()
{
	t_real x_mid = PSD_WIDTH/2 - 4;
	t_real y_mid = PSD_HEIGHT/2 + 4;
	t_real rad1 = 35.;
	t_real rad2 = 48.;
	t_real start_angle = M_PI*0.25;
	t_real end_angle = -M_PI*0.25;

	fs::path mask_file("mask.txt");
	if(fs::exists(mask_file))
	{
		std::ifstream mask(mask_file.string());
		if(mask)
		{
			mask
				>> x_mid >> y_mid
				>> rad1 >> rad2
				>> start_angle >> end_angle;

			start_angle = start_angle / 180. * M_PI;
			end_angle = end_angle / 180. * M_PI;
		}
	}

	bool ok = circular_mask("mask.png",
		x_mid, y_mid, rad1, rad2,
		start_angle, end_angle);

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

