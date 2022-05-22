/**
 * convert tof files to png images
 * @author Tobias Weber <tweber@ill.fr>
 * @date 15/nov/2021
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <bitset>

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


/**
 * extract image and count data from a tof file
 */
bool convert_tof(const fs::path& tof_file, const fs::path& out_file)
{
	// tof file data type for counts
	using t_data = std::uint32_t;

	// mask file
	using t_mask = std::bitset<PSD_HEIGHT*PSD_WIDTH>;
	fs::path mask_file("mask.png");
	gil::gray8_image_t mask;
	std::unique_ptr<t_mask> mask_bits;
	try
	{
		if(fs::exists(mask_file))
		{
			gil::read_image(mask_file.string(), mask, gil::png_tag());
			if(mask.width() == PSD_WIDTH && mask.height() == PSD_HEIGHT)
			{
				auto mask_view = gil::view(mask);
				mask_bits = std::make_unique<t_mask>();

				// load mask bits
				for(unsigned y=0; y<PSD_HEIGHT; ++y)
				{
					auto mask_row = mask_view.row_begin(y);

					for(unsigned x=0; x<PSD_WIDTH; ++x)
					{
						bool b = (*(mask_row + x) != 0);
						(*mask_bits)[y*PSD_WIDTH + x] = b;
					}
				}
			}
			else
			{
				std::cerr << "Error: Mask needs to be of size "
					<< PSD_WIDTH << " x " << PSD_HEIGHT
					<< "." << std::endl;
			}
		}
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Cannot load mask \""
			<< mask_file.string() << "\"."
			<< std::endl;
	}

	// total image, summing over all channels
	gil::gray16_image_t total_png(PSD_WIDTH, PSD_HEIGHT);
	auto total_view = gil::view(total_png);

	// output file for tof channel counts
	std::ostringstream cnts_file;
	cnts_file << out_file.string() << ".cnt";
	std::ofstream ofstr_cnts(cnts_file.str());
	ofstr_cnts << "#"
		<< std::setw(15) << std::right << "channel"
		<< std::setw(16) << std::right << "counts"
		<< std::setw(16) << std::right << "counts in mask"
		<< std::endl;

	// iterate over tof channels
	for(unsigned t=0; t<TOF_COUNT; ++t)
	{
		// output file for tof channel
		std::ostringstream ostr_file;
		ostr_file << out_file.string() << "_" << t << ".png";
		std::string png_file = ostr_file.str();

		ios::mapped_file_source file(tof_file,
			PSD_WIDTH*PSD_HEIGHT*sizeof(t_data),
			t*PSD_WIDTH*PSD_HEIGHT*sizeof(t_data));

		if(!file.is_open())
			return false;

		const t_data* data = reinterpret_cast<const t_data*>(file.data());

		// image of a tof channel
		gil::gray16_image_t png(PSD_WIDTH, PSD_HEIGHT);
		auto view = gil::view(png);

		std::uint64_t counts = 0;
		std::uint64_t counts_mask = 0;
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

				if(mask_bits)
				{
					bool in_mask = (*mask_bits)[y*PSD_WIDTH + x];
					if(in_mask)
						counts_mask += cnt;
				}
			}
		}

		file.close();

		ofstr_cnts
			<< std::setw(16) << std::right << t
			<< std::setw(16) << std::right << counts;
		if(mask_bits)
			ofstr_cnts << std::setw(16) << std::right << counts_mask;
		ofstr_cnts << std::endl;

		if(counts)
		{
			// write channel image
			gil::write_view(png_file, view, gil::png_tag());
		}
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
