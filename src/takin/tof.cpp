/**
 * convert tof files to png images
 * @author Tobias Weber <tweber@ill.fr>
 * @date 15/nov/2021
 * @license GPLv2 (see 'LICENSE' file)
 *
 * g++ -o tof tof.cpp -lboost_system -lboost_filesystem -lboost_iostreams -lpng
 */

#include <iostream>
#include <sstream>
#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
namespace ios = boost::iostreams;

#include <boost/gil/image.hpp>
#include <boost/gil/extension/io/png.hpp>
namespace gil = boost::gil;


#define PSD_WIDTH  128
#define PSD_HEIGHT 128
#define TOF_COUNT  128


bool convert_tof(const fs::path& tof_file, const fs::path& out_file)
{
	using t_data = std::uint32_t;

	gil::gray16_image_t total_png(PSD_WIDTH, PSD_HEIGHT);
	auto total_view = gil::view(total_png);

	for(unsigned t=0; t<TOF_COUNT; ++t)
	{
		std::ostringstream ostr_file;
		ostr_file << out_file.string() << "_" << t << ".png";
		std::string png_file = ostr_file.str();

		ios::mapped_file_source file(tof_file,
			PSD_WIDTH*PSD_HEIGHT*sizeof(t_data),
			t*PSD_WIDTH*PSD_HEIGHT*sizeof(t_data));

		if(!file.is_open())
			return false;

		const t_data* data = reinterpret_cast<const t_data*>(file.data());

		gil::gray16_image_t png(PSD_WIDTH, PSD_HEIGHT);
		auto view = gil::view(png);

		std::uint64_t counts = 0;
		for(unsigned y=0; y<PSD_HEIGHT; ++y)
		{
			auto row = view.row_begin(y);
			auto total_row = total_view.row_begin(y);

			for(unsigned x=0; x<PSD_WIDTH; ++x)
			{
				t_data cnt = data[y*PSD_WIDTH + x];

				*(row + x) = cnt;

				if(t == 0)
					*(total_row + x) = cnt;
				else
					*(total_row + x) = *(total_row + x) + cnt;

				counts += cnt;
				//std::cout << std::hex << cnt << " ";
			}
			//std::cout << std::endl;
		}

		file.close();

		if(counts)
		{
			//std::cout << "Channel " << t << ": " << std::dec << counts << std::endl;
			gil::write_view(png_file, view, gil::png_tag());
		}
	}

	std::string total_png_file = out_file.string() + ".png";
	gil::write_view(total_png_file, total_view, gil::png_tag());
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
			std::cerr << "File \"" << argv[i] << "\" does not exist!" << std::endl;
			continue;
		}


		fs::path file_out = file;
		file_out.replace_extension("");

		if(convert_tof(file, file_out))
			std::cout << "Converted file \"" << argv[i] << "\"." << std::endl;
		else
			std::cerr << "Failed to convert file \"" << argv[i] << "\"." << std::endl;
	}

	return 0;
}