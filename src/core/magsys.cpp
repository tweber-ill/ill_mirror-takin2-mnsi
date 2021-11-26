/**
 * magnetic system interface
 * @author Tobias Weber <tweber@ill.fr>
 * @date jul-18
 * @desc This file implements the theoretical skyrmion model by M. Garst and J. Waizner, references:
 *	- https://doi.org/10.1088/1361-6463/aa7573
 *	- https://kups.ub.uni-koeln.de/7937/
 *	- Personal communications with M. Garst, 2017-2020.
 * @desc This file is based on:
 *	- The descriptions and Mathematica implementations of the different skyrmion model versions by M. Garst and J. Waizner, 2016-2020,
 *	- The 2016 Python implementations by M. Kugler and G. Brandl of the first version of the skyrmion model.
 *	  This present version started as a C++ port of that Python implementation by M. Kugler and G. Brandl,
 *	  that was then adapted to new theoretical model revisions provided by M. Garst.
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "magsys.h"

#define USE_LAPACK
//#include "tlibs2/libs/math17.h"
#include "tlibs2-extras/math17.h"
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
	int iMaxOrder, bool bFixXR, bool bFixYR, bool bFixZR, bool bFixXI, bool bFixYI, bool bFixZI)
{
	auto& sys = *this;
	const auto& fourier = sys.GetFourier();

	std::vector<tl2::t_real_min> vars, errs;
	std::vector<std::string> names;
	std::vector<bool> fixed;
	int fourier_idx = 0;
	t_real err_def = 5.;

	for(const auto &vec : fourier)
	{
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

		++fourier_idx;
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
				{ args[i]+imag*args[i+1], args[i+2]+imag*args[i+3], args[i+4]+imag*args[i+5] });
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
	bool bFixXR, bool bFixYR, bool bFixZR, bool bFixXI, bool bFixYI, bool bFixZI) const
{
	bool calc_angles = 1;

	std::ofstream ofstr(file);
	std::shared_ptr<std::ofstream> ofstrGpl, ofstrGplBc2;

	if(calc_angles)
	{
		ofstrGpl.reset(new std::ofstream(std::string(file)+"_ang.gpl"));
		ofstrGplBc2.reset(new std::ofstream(std::string(file)+"_Bc2.gpl"));

		ofstrGpl->precision(8);
		ofstrGplBc2->precision(8);

		(*ofstrGpl) << "#!/usr/local/bin/gnuplot -p\n\n";
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
	for(t_real B = 0; B < 200; B += 2.5)
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


	std::vector<std::vector<t_real>> Bc2s;
	Bc2s.resize(Ts.size());

	std::size_t iState = 0;
	for(t_real B : Bs)
	{
		if(ofstrGpl) (*ofstrGpl) << std::setw(16) << B << " ";

		for(std::size_t T_idx=0; T_idx<Ts.size(); ++T_idx)
		{
			t_real T = Ts[T_idx];

			mag->SetT(T);
			mag->SetB(B);

			bool ok = mag->minimise(iMaxOrder, bFixXR, bFixYR, bFixZR, bFixXI, bFixYI, bFixZI);
			if(T_idx == 0)	// run twice to get better initial values
				ok = mag->minimise(iMaxOrder, bFixXR, bFixYR, bFixZR, bFixXI, bFixYI, bFixZI);
			const auto& fourier = mag->GetFourier();
			auto f = mag->F();

			const std::string labState = "state_" + tl2::var_to_str(iState);
			ofstr << "\t<" << labState << ">\n";

			t_real ang = std::atan2(std::norm(fourier[1][0])*std::sqrt(t_real(2)), std::norm(fourier[0][2]));
			if(ok && ang < tl2::pi<t_real>/4. && ang > 0.01)
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

			std::cout << "\rB = " << B << ", T = " << T << "        ";
			std::cout.flush();

			ofstr << "\t</" << labState << ">\n";
			++iState;
		}

		if(ofstrGpl) (*ofstrGpl) << "\n";
	}

	if(ofstrGpl) (*ofstrGpl) << "e\n";
	std::cout << std::endl;
	ofstr << "</states>\n";


	if(ofstrGplBc2)
	{
		std::ostringstream ostrData;
		ostrData.precision(ofstrGplBc2->precision());
		for(std::size_t T_idx=0; T_idx<Ts.size(); ++T_idx)
		{
			t_real T = Ts[T_idx];

			t_real mean = tl2::mean_value(Bc2s[T_idx]);
			t_real dev = tl2::std_dev(Bc2s[T_idx]);

			ostrData << std::setw(16) << T << " ";
			ostrData << std::setw(16) << mean << " ";
			ostrData << std::setw(16) << dev << "\n";
		}
		ostrData << "e\n";

		// fit
		(*ofstrGplBc2) << "amp = 1\nex = 0.5\n";
		(*ofstrGplBc2) << "func(T) = amp * (-T)**ex\n\n";

		(*ofstrGplBc2) << "fit func(x) \"-\" u 1:2:3 yerr via amp, ex\n";
		(*ofstrGplBc2) << ostrData.str() << "\n";

		// plot
		(*ofstrGplBc2) << "set xlabel \"T\"\nset ylabel \"Bc2\"\n";
		(*ofstrGplBc2) << "set autosc xfix\nset autosc yfix\n\n";

		(*ofstrGplBc2) << "plot \"-\" u 1:2:3 w yerrorbars pt 7, func(x) w lines lw 2\n";
		(*ofstrGplBc2) << ostrData.str();
	}

	return true;
}


#else


template<class t_real, class t_cplx, int ORDER_FOURIER>
bool MagSystem<t_real, t_cplx, ORDER_FOURIER>::minimise(
	int iMaxOrder, bool bFixXR, bool bFixYR, bool bFixZR, bool bFixXI, bool bFixYI, bool bFixZI)
{
	return 0;
}


template<class t_real, class t_cplx, int ORDER_FOURIER>
bool MagSystem<t_real, t_cplx, ORDER_FOURIER>::SaveStates(
	const char *file, int iMaxOrder,
	bool bFixXR, bool bFixYR, bool bFixZR, bool bFixXI, bool bFixYI, bool bFixZI) const
{
	return 0;
}


#endif
