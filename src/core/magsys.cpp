/**
 * magnetic system interface
 * @author Tobias Weber <tweber@ill.fr>
 * @date jul-18
 * @desc This file implements the theoretical skyrmion model by M. Garst and J. Waizner, references:
 *      - M. Garst, J. Waizner, and D. Grundler, J. Phys. D: Appl. Phys. 50, 293002 (2017), https://doi.org/10.1088/1361-6463/aa7573
 *      - J. Waizner, PhD thesis (2017), Universität zu Köln, https://kups.ub.uni-koeln.de/7937/
 *	- Personal communications with M. Garst, 2017-2020.
 * @desc This file is based on:
 *	- The descriptions and Mathematica implementations of the different skyrmion model versions by M. Garst and J. Waizner, 2016-2020,
 *	- The 2016 Python implementations by M. Kugler and G. Brandl of the first version of the skyrmion model.
 *	  This present version started as a C++ port of M. Kugler's and G. Brandl's Python implementation,
 *	  that was then adapted to new theoretical model revisions provided by M. Garst.
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "magsys.h"

#define USE_LAPACK
#include "tlibs2/libs/math17.h"
#include "tlibs2/libs/str.h"
#include "tlibs2/libs/log.h"

#ifndef NO_MINIMISATION
	#include "tlibs2/libs/fit.h"
#endif

#include <fstream>


// instantiation
#ifdef DEF_SKX_ORDER
	template class MagSystem<double, std::complex<double>, (DEF_SKX_ORDER+1)*DEF_SKX_ORDER/2>;
#endif
#ifdef DEF_HELI_ORDER
	template class MagSystem<double, std::complex<double>, DEF_HELI_ORDER>;
#endif


#ifndef NO_MINIMISATION


/**
 * minimisation of free energy
 */
template<class t_real, class t_cplx, int ORDER_FOURIER>
bool MagSystem<t_real, t_cplx, ORDER_FOURIER>::minimise(
	int iMaxOrder,
	bool bFixXR, bool bFixYR, bool bFixZR,
	bool bFixXI, bool bFixYI, bool bFixZI)
{
	auto& sys = *this;
	const auto& fourier = sys.GetFourier();

	std::vector<tl2::t_real_min> vars, errs;
	std::vector<std::string> names;
	std::vector<bool> fixed;
	t_real err_def = 5.;

	for(int fourier_idx=0; fourier_idx<ORDER_FOURIER+1; ++fourier_idx)
	{
		const auto &vec = fourier[fourier_idx];

		vars.push_back(vec[0].real()); vars.push_back(vec[0].imag());
		vars.push_back(vec[1].real()); vars.push_back(vec[1].imag());
		vars.push_back(vec[2].real()); vars.push_back(vec[2].imag());

		errs.push_back(err_def); errs.push_back(err_def);
		errs.push_back(err_def); errs.push_back(err_def);
		errs.push_back(err_def); errs.push_back(err_def);

		names.push_back(std::string("m") + tl2::var_to_str(fourier_idx) + std::string("_x_real"));
		names.push_back(std::string("m") + tl2::var_to_str(fourier_idx) + std::string("_x_imag"));
		names.push_back(std::string("m") + tl2::var_to_str(fourier_idx) + std::string("_y_real"));
		names.push_back(std::string("m") + tl2::var_to_str(fourier_idx) + std::string("_y_imag"));
		names.push_back(std::string("m") + tl2::var_to_str(fourier_idx) + std::string("_z_real"));
		names.push_back(std::string("m") + tl2::var_to_str(fourier_idx) + std::string("_z_imag"));

		/*if(fourier_idx > iMaxOrder)
		{
			fixed.push_back(1); fixed.push_back(1);
			fixed.push_back(1); fixed.push_back(1);
			fixed.push_back(1); fixed.push_back(1);
		}
		else*/ if(fourier_idx == 0)
		{
			fixed.push_back(1); fixed.push_back(1);
			fixed.push_back(1); fixed.push_back(1);
			fixed.push_back(bFixZR); fixed.push_back(1);
		}
		else
		{
			fixed.push_back(bFixXR); fixed.push_back(bFixXI);
			fixed.push_back(bFixYR); fixed.push_back(bFixYI);
			fixed.push_back(bFixZR); fixed.push_back(bFixZI);
		}
	}


	// convert minimiser args to complex vector
	auto vars_to_fourier = [](const std::vector<tl2::t_real_min>& args)
		-> std::vector<ublas::vector<t_cplx>>
	{
		constexpr auto imag = t_cplx(0,1);

		std::vector<ublas::vector<t_cplx>> _fourier;
		for(std::size_t i=0; i<args.size(); i+=6)
		{
			auto m = tl2::make_vec<ublas::vector<t_cplx>>(
			{
				args[i] + imag*args[i+1],
				args[i+2] + imag*args[i+3],
				args[i+4] + imag*args[i+5]
			});
			_fourier.emplace_back(std::move(m));
		}
		return _fourier;
	};


	// lambda function to minimise
	auto func = [this, &vars_to_fourier](const auto& ...pack) -> tl2::t_real_min
	{
		const std::vector<tl2::t_real_min>& args = {{ pack... }};

		auto newsys = this->copyCastSys();
		newsys->SetFourier(vars_to_fourier(args));

		return newsys->F();
	};

	constexpr std::size_t NUM_MINPARAMS = 6*(ORDER_FOURIER+1);
	if(vars.size() != NUM_MINPARAMS)
	{
		tl2::log_err("Fourier component number mismatch: ", vars.size(), " != ", NUM_MINPARAMS, ".");
		return false;
	}

	bool ok = tl2::minimise<tl2::t_real_min, NUM_MINPARAMS>(func, names, vars, errs, &fixed, m_debug);
	sys.SetFourier(vars_to_fourier(vars));
	return ok;
}



/**
 * calculates and saves states for phase diagram
 */
template<class t_real, class t_cplx, int ORDER_FOURIER>
bool MagSystem<t_real, t_cplx, ORDER_FOURIER>::SaveStates(
	const char *file, int iMaxOrder,
	bool bFixXR, bool bFixYR, bool bFixZR,
	bool bFixXI, bool bFixYI, bool bFixZI) const
{
	bool calc_angles = 1;

	std::ofstream ofstr(file);
	std::shared_ptr<std::ofstream> ofstrGpl, ofstrGplF, ofstrGplBc2;

	if(calc_angles)
	{
		ofstrGpl.reset(new std::ofstream(std::string(file)+"_ang.gpl"));
		ofstrGplF.reset(new std::ofstream(std::string(file)+"_F.gpl"));
		ofstrGplBc2.reset(new std::ofstream(std::string(file)+"_Bc2.gpl"));

		ofstrGpl->precision(8);
		ofstrGplF->precision(8);
		ofstrGplBc2->precision(8);

		(*ofstrGpl) << "#!/usr/local/bin/gnuplot -p\n\n";
		(*ofstrGplF) << "#!/usr/local/bin/gnuplot -p\n\n";
		(*ofstrGplBc2) << "#!/usr/local/bin/gnuplot -p\n\n";
	}


	if(!ofstr)
		return false;
	ofstr.precision(8);

	ofstr << "<states>\n";
	ofstr << "\t<order> " << ORDER_FOURIER << " </order>\n";
	ofstr << "\t<chi> " << g_chi<t_real> << " </chi>\n";
	ofstr << "\t<hoc> " << g_hoc<t_real> << " </hoc>\n";


	auto mag = this->copyCastSys();

	std::vector<t_real> Ts, Bs;
	for(t_real T = -20000; T < 0; T += 250.)
		Ts.push_back(T);
	for(t_real B = 0; B < 200; B += 1.0)
		Bs.push_back(B);


	if(ofstrGpl)
	{
		(*ofstrGpl) << "set xlabel \"T\"\nset ylabel \"B\"\n";
		(*ofstrGpl) << "set autosc xfix\nset autosc yfix\nset cbrange [0 : pi/2]\n\n";
		(*ofstrGpl) << "plot \"-\" mat nonuni w ima\n";
		(*ofstrGpl) << std::setw(16) << "0" << " ";
		for(t_real T : Ts)
			(*ofstrGpl) << std::setw(16) << T << " ";
		(*ofstrGpl) << "\n";
	}

	if(ofstrGplF)
	{
		(*ofstrGplF) << "set xlabel \"T\"\nset ylabel \"B\"\n";
		(*ofstrGplF) << "set autosc xfix\nset autosc yfix\n\n";
		(*ofstrGplF) << "plot \"-\" mat nonuni w ima\n";
		(*ofstrGplF) << std::setw(16) << "0" << " ";
		for(t_real T : Ts)
			(*ofstrGplF) << std::setw(16) << T << " ";
		(*ofstrGplF) << "\n";
	}


	std::vector<std::vector<t_real>> Bc2s;
	Bc2s.resize(Ts.size());

	std::size_t iState = 0;
	for(t_real B : Bs)
	{
		if(ofstrGpl) (*ofstrGpl) << std::setw(16) << B << " ";
		if(ofstrGplF) (*ofstrGplF) << std::setw(16) << B << " ";

		for(std::size_t T_idx=0; T_idx<Ts.size(); ++T_idx)
		{
			t_real T = Ts[T_idx];

			mag->SetT(T, false);
			mag->SetB(B, false);

			bool ok = mag->minimise(iMaxOrder, bFixXR, bFixYR, bFixZR, bFixXI, bFixYI, bFixZI);
			if(T_idx == 0)	// run twice to get better initial values
				ok = mag->minimise(iMaxOrder, bFixXR, bFixYR, bFixZR, bFixXI, bFixYI, bFixZI);
			const auto& fourier = mag->GetFourier();
			t_real f = mag->F();

			const std::string labState = "state_" + tl2::var_to_str(iState);
			ofstr << "\t<" << labState << ">\n";

			const t_real ang_min = 0.01;
			t_real ang = std::atan2(std::abs(2.*fourier[1][0]),
				std::abs(fourier[0][2]));
			if(ok && tl2::is_in_angular_range<t_real>(ang_min, tl2::pi<t_real>/4.-ang_min, ang))
			{
				t_real Bc2 = B / std::cos(ang);
				Bc2s[T_idx].push_back(Bc2);

				ofstr << "\t\t<ang> " << ang << " </ang>\n";
				ofstr << "\t\t<Bc2> " << Bc2 << " </Bc2>\n";
			}

			ofstr << "\t\t<ok> " << int(ok) << " </ok>\n";
			ofstr << "\t\t<B> " << B << " </B>\n";
			ofstr << "\t\t<T> " << T << " </T>\n";
			ofstr << "\t\t<F> " << f << " </F>\n";
			for(std::size_t i=0; i<fourier.size(); ++i)
			{
				const std::string labM = "M_" + tl2::var_to_str(i);
				ofstr << "\t\t<" << labM << "> "
					<< fourier[i][0] << ";\t" << fourier[i][1] << ";\t" << fourier[i][2]
					<< " </" << labM << ">\n";
			}


			if(ofstrGpl) (*ofstrGpl) << std::setw(16) << ang << " ";
			if(ofstrGplF) (*ofstrGplF) << std::setw(16) << f << " ";

			std::cout << "\rB = " << B << ", T = " << T << "        ";
			std::cout.flush();

			ofstr << "\t</" << labState << ">\n";
			++iState;
		}

		if(ofstrGpl) (*ofstrGpl) << "\n";
		if(ofstrGplF) (*ofstrGplF) << "\n";
	}

	if(ofstrGpl) (*ofstrGpl) << "e\n";
	if(ofstrGplF) (*ofstrGplF) << "e\n";
	std::cout << std::endl;
	ofstr << "</states>\n";


	if(ofstrGplBc2)
	{
		// data
		std::ostringstream ostrData;
		ostrData.precision(ofstrGplBc2->precision());
		ostrData << "$data << ENDDATA\n";
		for(std::size_t T_idx=0; T_idx<Ts.size(); ++T_idx)
		{
			t_real T = Ts[T_idx];

			t_real mean = tl2::mean_value(Bc2s[T_idx]);
			t_real dev = tl2::std_dev(Bc2s[T_idx]);
			if(tl2::float_equal(dev, t_real(0)))
				continue;

			ostrData << std::setw(16) << T << " ";
			ostrData << std::setw(16) << mean << " ";
			ostrData << std::setw(16) << dev << "\n";
		}
		ostrData << "ENDDATA\n";
		(*ofstrGplBc2) << ostrData.str() << "\n";

		// fit
		(*ofstrGplBc2) << "func1(T, amp, ex) = amp * (-T)**ex\n";
		(*ofstrGplBc2) << "func2(T, amp) = 2 * amp * sqrt(-0.5 - 0.5*T)\n\n";
		(*ofstrGplBc2) << "amp1 = 1\nex1 = 0.5\n";
		(*ofstrGplBc2) << "amp2 = 1\n\n";
		(*ofstrGplBc2) << "fit func1(x, amp1, ex1) \"$data\" u 1:2:3 yerr via amp1, ex1\n";
		(*ofstrGplBc2) << "fit func2(x, amp2) \"$data\" u 1:2:3 yerr via amp2\n\n";

		// plot
		(*ofstrGplBc2) << "set xlabel \"T\"\nset ylabel \"Bc2\"\n";
		(*ofstrGplBc2) << "set autosc xfix\nset autosc yfix\n\n";
		(*ofstrGplBc2) << "plot \"$data\" u 1:2:3 w yerrorbars pt 7, "
			<< "func1(x, amp1, ex1) w lines lw 2, "
			<< "func2(x, amp2) lw 2\n\n";

		(*ofstrGplBc2) << "msg = sprintf(\"\\nModel 1:\\n\\tamp = %.8f\\n\\tex = %.8f\\n\\nModel 2:\\n\\tamp = %.8f\\n\", amp1, ex1, amp2)\n";
		(*ofstrGplBc2) << "print(msg)\n";
	}

	return true;
}


#else  // use dummy implementations if not needed


template<class t_real, class t_cplx, int ORDER_FOURIER>
bool MagSystem<t_real, t_cplx, ORDER_FOURIER>::minimise(
	int /*iMaxOrder*/,
	bool /*bFixXR*/, bool /*bFixYR*/, bool /*bFixZR*/,
	bool /*bFixXI*/, bool /*bFixYI*/, bool /*bFixZI*/)
{
	return 0;
}


template<class t_real, class t_cplx, int ORDER_FOURIER>
bool MagSystem<t_real, t_cplx, ORDER_FOURIER>::SaveStates(
	const char* /*file*/, int /*iMaxOrder*/,
	bool /*bFixXR*/, bool /*bFixYR*/, bool /*bFixZR*/,
	bool /*bFixXI*/, bool /*bFixYI*/, bool /*bFixZI*/) const
{
	return 0;
}


#endif
