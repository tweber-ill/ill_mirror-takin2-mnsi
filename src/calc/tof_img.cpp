/**
 * converts tof files to png images
 * @author Tobias Weber <tweber@ill.fr>
 * @date 15/nov/2021
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <bitset>
#include <tuple>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
namespace ios = boost::iostreams;

#include <boost/gil/image.hpp>
#include <boost/gil/extension/io/png.hpp>
namespace gil = boost::gil;

#include "tof.h"

//#define WRITE_TOF_CHANNEL


/**
 * extract image data from a tof file
 */
bool convert_tof(const fs::path& tof_file, const fs::path& out_file)
{
	// tof file data type for counts
	using t_data = std::uint32_t;

	// total image, summing over all channels
	gil::gray16_image_t total_png(PSD_WIDTH, PSD_HEIGHT);
	auto total_view = gil::view(total_png);

	// iterate over tof channels
	for(unsigned t=0; t<TOF_COUNT; ++t)
	{
		ios::mapped_file_source file(tof_file,
			PSD_WIDTH*PSD_HEIGHT*sizeof(t_data),
			t*PSD_WIDTH*PSD_HEIGHT*sizeof(t_data));
		if(!file.is_open())
			return false;

		const t_data* data = reinterpret_cast<const t_data*>(file.data());

#ifdef WRITE_TOF_CHANNEL
		// output file for tof channel
		std::ostringstream ostr_file;
		ostr_file << out_file.string() << "_" << t << ".png";
		std::string png_file = ostr_file.str();

		// image of a tof channel
		gil::gray16_image_t png(PSD_WIDTH, PSD_HEIGHT);
		auto view = gil::view(png);
#endif

		std::uint64_t counts = 0;
		for(unsigned y=0; y<PSD_HEIGHT; ++y)
		{
#ifdef WRITE_TOF_CHANNEL
			auto row = view.row_begin(y);
#endif
			auto total_row = total_view.row_begin(y);

			for(unsigned x=0; x<PSD_WIDTH; ++x)
			{
				t_data cnt = data[y*PSD_WIDTH + x];
#ifdef WRITE_TOF_CHANNEL
				*(row + x) = cnt;
#endif

				if(t == 0)
					*(total_row + x) = cnt;
				else
					*(total_row + x) = *(total_row + x) + cnt;

				counts += cnt;
			}
		}

		file.close();

#ifdef WRITE_TOF_CHANNEL
		if(counts)
		{
			// write channel image
			gil::write_view(png_file, view, gil::png_tag());
		}
#endif
	}

	// write total image
	std::string total_png_file = out_file.string() + ".png";
	gil::write_view(total_png_file, total_view, gil::png_tag());

	// meta information at the end of the file
	std::size_t size = fs::file_size(tof_file);
	std::ptrdiff_t size_rest = size - PSD_WIDTH*PSD_HEIGHT*TOF_COUNT*sizeof(t_data);
	if(size_rest > 0)
	{
		ios::mapped_file_source file(tof_file,
			size_rest,
			PSD_WIDTH*PSD_HEIGHT*TOF_COUNT*sizeof(t_data));

		if(file.is_open())
		{
			const char* data = reinterpret_cast<const char*>(file.data());

			std::string txt_file = out_file.string() + ".txt";
			std::ofstream ofstr_txt(txt_file);

			ofstr_txt.write(data, size_rest);
		}
	}

	return true;
}


int main(int argc, char** argv)
{
	if(argc < 2)
	{
		std::cerr << "Please give some TOF files." << std::endl;
		return -1;
	}

	for(int i=1; i<argc; ++i)
	{
		fs::path file(argv[i]);

		if(!fs::exists(file))
		{
			std::cerr << "Error: File \"" << file.string()
				<< "\" does not exist!" << std::endl;
			continue;
		}


		fs::path file_out = file.filename();
		file_out.replace_extension("");

		if(convert_tof(file, file_out))
		{
			std::cerr << "Converted file \""
				<< file.string() << "\"." << std::endl;
		}
		else
		{
			std::cerr << "Error: Failed to process file \""
				<< file.string() << "\"." << std::endl;
		}
	}

	return 0;
}
