/**
 * calculates the polarisation of tof files
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

#include "tlibs2/libs/fit.h"
#include "tof.h"

using t_real = tl2::t_real_min;


/**
 * extract image and count data from a tof file
 */
std::tuple<bool, t_real, t_real, t_real>
process_tof(const fs::path& tof_file, const fs::path& out_file, int foil=3)
{
	// tof file data type for counts
	using t_data = std::uint32_t;

	// image mask file giving allowed pixels as non-black value
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
		std::cerr << "Error: Cannot load mask \""
			<< mask_file.string() << "\"."
			<< std::endl;
	}

	// tof mask file giving tof foil start and end indices
	unsigned tof_mask_start = foil * CHANNELS_PER_FOIL;
	unsigned tof_mask_end = (foil+1) * CHANNELS_PER_FOIL;

	// output file for tof channel counts
	std::ostringstream cnts_file;
	cnts_file << out_file.string() << ".cnt";
	std::ofstream ofstr_cnts(cnts_file.str());
	ofstr_cnts << "#"
		<< std::setw(15) << std::right << "channel"
		<< std::setw(16) << std::right << "counts"
		<< std::setw(16) << std::right << "counts in mask"
		<< std::endl;

	unsigned tof_channel = 0;
	std::vector<t_real> tof_channels;
	std::vector<t_real> tof_counts;
	std::vector<t_real> tof_errors;

	// iterate over tof channels
	for(unsigned t=0; t<TOF_COUNT; ++t)
	{
		ios::mapped_file_source file(tof_file,
			PSD_WIDTH*PSD_HEIGHT*sizeof(t_data),
			t*PSD_WIDTH*PSD_HEIGHT*sizeof(t_data));
		if(!file.is_open())
			return std::make_tuple(false, -1., -1., -1.);

		const t_data* data = reinterpret_cast<const t_data*>(file.data());

		std::uint64_t counts = 0;
		std::uint64_t counts_mask = 0;
		for(unsigned y=0; y<PSD_HEIGHT; ++y)
		{
			for(unsigned x=0; x<PSD_WIDTH; ++x)
			{
				t_data cnt = data[y*PSD_WIDTH + x];
				counts += cnt;

				// counts inside mask
				bool in_mask = true;
				if(mask_bits)
					in_mask = (*mask_bits)[y*PSD_WIDTH + x];
				if(in_mask && t>=tof_mask_start && t<tof_mask_end)
					counts_mask += cnt;
			}
		}

		// counts inside tof mask
		if(t>=tof_mask_start && t<tof_mask_end)
		{
			tof_channels.push_back(tof_channel++);
			tof_counts.push_back(counts_mask);
			tof_errors.push_back(counts_mask ? std::sqrt(counts_mask) : 1);
		}

		file.close();

		ofstr_cnts
			<< std::setw(16) << std::right << t
			<< std::setw(16) << std::right << counts
			<< std::setw(16) << std::right << counts_mask;
		ofstr_cnts << std::endl;
	}

	// get spin-echo time from meta information at the end of the file
	std::size_t size = fs::file_size(tof_file);
	std::ptrdiff_t size_rest = size - PSD_WIDTH*PSD_HEIGHT*TOF_COUNT*sizeof(t_data);
	t_real tau = 0.;  // spin-echo time in ps
	if(size_rest > 0)
	{
		ios::mapped_file_source file(tof_file,
			size_rest,
			PSD_WIDTH*PSD_HEIGHT*TOF_COUNT*sizeof(t_data));

		if(file.is_open())
		{
			const char* data = reinterpret_cast<const char*>(file.data());

			//std::string txt_file = out_file.string() + ".txt";
			//std::ofstream ofstr_txt(txt_file);
			//ofstr_txt.write(data, size_rest);

			// get spin-echo time from metadata
			std::string tau_key = "echotime_value";
			const char* tau_ptr = std::search(data, data+size_rest, tau_key.begin(), tau_key.end());
			if(tau_ptr != data+size_rest)
			{
				// sometimes the key exists multiple times with the second one being correct
				const char* tau_ptr_2 = std::search(tau_ptr+1, data+size_rest, tau_key.begin(), tau_key.end());
				if(tau_ptr_2 != data+size_rest)
					tau_ptr = tau_ptr_2;

				const char *tau_ptr_end = std::find(tau_ptr, data+size_rest, '\n');
				std::string tau_str(tau_ptr, tau_ptr_end);
				std::istringstream istr_tau(tau_str);
				std::string tau_key, tau_dummy, tau_unit;
				istr_tau >> tau_key >> tau_dummy >> tau >> tau_unit;
				if(tau_unit == "ns")
					tau *= 1000.;
				else if(tau_unit == "us")
					tau *= 1e6;
				else if(tau_unit == "ms")
					tau *= 1e9;
				else if(tau_unit == "s")
					tau *= 1e12;
			}
		}
	}

	// fit spin-echo sine
	auto se_sine = [](t_real x, t_real amp, t_real freq,
		t_real phase, t_real offs)
	{
		return amp*std::sin(freq*x + phase) + offs;
	};

	// remove invalid first and last point
	tof_channels.erase(tof_channels.begin()); tof_channels.pop_back();
	tof_counts.erase(tof_counts.begin()); tof_counts.pop_back();
	tof_errors.erase(tof_errors.begin()); tof_errors.pop_back();

	t_real tof_counts_mean = tl2::mean_value(tof_counts);
	t_real tof_counts_dev = tl2::std_dev(tof_counts);
	//std::cout << "tof counts: " << tof_counts_mean << " +- " << tof_counts_dev << std::endl;

	std::vector<std::string> fit_vars{{"amp", "freq", "phase", "offs"}};
	std::vector<bool> fixed_vars{{ false, true, false, false }};
	std::vector<t_real> fit_vals{{
		tof_counts_dev,
		t_real(2.)*tl2::pi<t_real> / t_real(/*tof_channels.size()*/ CHANNELS_PER_FOIL),
		0.,
		tof_counts_mean
	}};
	std::vector<t_real> fit_errs{{
		tof_counts_dev / t_real(5.),
		t_real(10. / 180.) * tl2::pi<t_real>,
		tl2::pi<t_real>,
		tof_counts_mean / t_real(5.)
	}};
	bool fit_ok = tl2::fit<t_real, 5>(se_sine,
		tof_channels, tof_counts, tof_errors,
		fit_vars, fit_vals, fit_errs, &fixed_vars, false);

	t_real pol = fit_vals[0] / fit_vals[3];
	t_real pol_err = std::sqrt(std::pow(fit_errs[0]/fit_vals[3], 2.)
		+ std::pow(fit_vals[0]*fit_errs[3]/(fit_vals[3]*fit_vals[3]), 2.));

	ofstr_cnts << "#\n"
		<< "# se sine fit valid: " << std::boolalpha << fit_ok << "\n"
		<< "# se sine fit: " << fit_vals[0] << " * sin(" << fit_vals[1]
			<< "*x + " << fit_vals[2] << ") + " << fit_vals[3] << "\n"
		<< "# se polarisation: " << pol << "\n"
		<< "# se polarisation error: " << pol_err  << "\n"
		<< "# plot with: gnuplot -p -e \"plot \\\"" << cnts_file.str()
			<< "\\\" u 1:3:(sqrt(\\$3)) w yerrorbars pt 7, " << fit_vals[0]
			<< " * sin(" << fit_vals[1] << "*x + " << fit_vals[2] << ") + "
			<< fit_vals[3] << "\"\n"
		<< "#" << std::endl;

	return std::make_tuple(fit_ok, tau, pol, pol_err);
}


int main(int argc, char** argv)
{
	if(argc < 2)
	{
		std::cerr << "Please give some TOF files." << std::endl;
		return -1;
	}

	std::cout.precision(6);
	std::cout
		<< std::setw(1) << std::right << "#"
		<< std::setw(4) << std::right << "idx"
		<< std::setw(5) << std::right << "foil"
		<< std::setw(48) << std::right << "file"
		<< std::setw(12) << std::right << "tau (ps)"
		<< std::setw(12) << std::right << "pol."
		<< std::setw(12) << std::right << "err."
		<< std::endl;

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

		for(int foil=0; foil<NUM_FOILS; ++foil)
		{
			bool ok = false;
			t_real tau = 0.;
			t_real pol = 0., pol_err = -1.;
			std::tie(ok, tau, pol, pol_err) = process_tof(file, file_out, foil);
			if(!ok)
			{
				std::cerr << "Error: Failed to process file \""
					<< file << "\", foil " << foil
					<< "." << std::endl;

			}

			std::ostringstream foil_data;
			foil_data << "pol_" << (foil+1) << ".dat";
			std::ofstream ofstr(foil_data.str(), std::ios_base::app);
			std::ostream* ostrs[]{ &std::cout, &ofstr };
			for(std::ostream* ostr : ostrs)
			{
				(*ostr)
					<< std::setw(5) << std::right << i
					<< std::setw(5) << std::right << (foil+1)
					<< std::setw(48) << std::right << file.string()
					<< std::setw(12) << std::right << tau
					<< std::setw(12) << std::right << pol
					<< std::setw(12) << std::right << pol_err
					<< std::endl;
			}
		}
	}

	return 0;
}
