/**
 * unites tof files
 * @author Tobias Weber <tweber@ill.fr>
 * @date 15/nov/2021
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <iostream>
#include <fstream>
#include <string>
#include <array>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
namespace ios = boost::iostreams;


#define PSD_WIDTH  128
#define PSD_HEIGHT 128
#define TOF_COUNT  128


/**
 * extract image and count data from a tof file
 */
bool unite_tof(const fs::path& united, const std::vector<fs::path>& tofs)
{
	// tof file data type for counts
	using t_data = std::uint32_t;

	std::ofstream ofstr(united.string());
	if(!ofstr)
		return false;
	auto sum_data = std::make_unique<std::array<t_data, TOF_COUNT*PSD_WIDTH*PSD_HEIGHT>>();
	std::memset(sum_data->data(), 0, TOF_COUNT*PSD_WIDTH*PSD_HEIGHT*sizeof(t_data));

	// unite count data
	for(const fs::path& tof_file : tofs)
	{
		ios::mapped_file_source file(tof_file,
			TOF_COUNT*PSD_WIDTH*PSD_HEIGHT*sizeof(t_data), 0);

		if(!file.is_open())
			return false;

		const t_data* data = reinterpret_cast<const t_data*>(file.data());
		for(std::size_t i=0; i<TOF_COUNT*PSD_WIDTH*PSD_HEIGHT; ++i)
			(*sum_data)[i] += data[i];

		file.close();
	}

	// write count data
	ofstr.write(reinterpret_cast<char*>(sum_data->data()),
		TOF_COUNT*PSD_WIDTH*PSD_HEIGHT*sizeof(t_data));


	// add the meta information of one of the input files
	const fs::path& tof_file = *tofs.begin();
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
			ofstr.write(data, size_rest);
		}
	}

	return true;
}


int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << argv[0]
			<< " <united.tof> <in1.tof> [in2.tof] [...]"
			<< std::endl;
		return -1;
	}

	fs::path united(argv[1]);
	std::vector<fs::path> files;

	for(int i=2; i<argc; ++i)
	{
		fs::path file(argv[i]);

		if(!fs::exists(file))
		{
			std::cerr << "Error: File \"" << file.string() << "\" does not exist!" << std::endl;
			continue;
		}

		files.emplace_back(std::move(file));
	}

	if(unite_tof(united, files))
	{
		std::cout << "Created file \"" << united.string() << "\"." << std::endl;
	}
	else
	{
		std::cerr << "Error creating file \"" << united.string() << "\"." << std::endl;
		return -1;
	}

	return 0;
}
