/**
 * Takin module
 * @author Tobias Weber <tweber@ill.fr>
 * @date sep-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#ifdef PLUGIN_APPLI
	#include "takin/tools/monteconvo/sqw_proc.h"
	#include "takin/tools/monteconvo/sqw_proc_impl.h"
#else
	#include <boost/dll/alias.hpp>
#endif

#include "takin/tools/monteconvo/sqwbase.h"
#include "takin/libs/version.h"

#include "core/skx.h"
#include "core/fp.h"
#include "core/heli.h"
#include "core/longfluct.h"

#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/str.h"


#ifndef DEF_SKX_ORDER
	#define SKX_ORDER 4
#else
	#define SKX_ORDER DEF_SKX_ORDER
#endif

#ifndef DEF_HELI_ORDER
	#define HELI_ORDER 4
#else
	#define HELI_ORDER DEF_HELI_ORDER
#endif


class SqwMod : public SqwBase
{
public:
	using SqwBase::t_var;
	using t_real = t_real_reso;
	using t_vec = ublas::vector<t_real>;
	using t_cplx = std::complex<t_real>;

protected:
	Skx<t_real, t_cplx, SKX_ORDER> m_skx;
	FP<t_real, t_cplx> m_fp;
	Heli<t_real, t_cplx, HELI_ORDER> m_heli;
	Longfluct m_lf;

	t_real m_dErange = 2.;	// E range around queried E to take into account
	t_real m_dSigma = t_real(0.05);
	t_real m_dS0 = t_real(1.);
	t_real m_dIncAmp = t_real(0.);
	t_real m_dIncSigma = t_real(0.05);
	t_real m_dT = 29.;
	t_real m_dB = 0.35;
	t_real m_dcut = 0.02;

	int m_iOnlyMode = -1;
	int m_iPolChan = 0;
	int m_iwhich_disp = 0;	// 0: skx, 1: fp, 2: heli
	int m_iProjNeutron = 1;
	int m_ionlylf = 0;

	t_vec m_vecG = tl2::make_vec<t_vec>({1,1,0});
	t_vec m_vecB = tl2::make_vec<t_vec>({1,1,0});
	t_vec m_vecPin = tl2::make_vec<t_vec>({1,-1,0});

public:
	SqwMod();
	SqwMod(const std::string& strCfgFile);
	virtual ~SqwMod();

	virtual std::tuple<std::vector<t_real>, std::vector<t_real>> disp(t_real dh, t_real dk, t_real dl) const override;
	virtual t_real operator()(t_real dh, t_real dk, t_real dl, t_real dE) const override;

	virtual std::vector<t_var> GetVars() const override;
	virtual void SetVars(const std::vector<t_var>&) override;
	virtual bool SetVarIfAvail(const std::string& strKey, const std::string& strNewVal) override;

	virtual SqwBase* shallow_copy() const override;
};


using t_real = typename SqwMod::t_real;
#include "core/skx_default_gs.cxx"


// ----------------------------------------------------------------------------
SqwMod::SqwMod()
{
	tl2::log_info("--------------------------------------------------------------------------------");
	tl2::log_info("This is the Takin MnSi magnon dynamics module,");
	tl2::log_info("\tby T. Weber <tweber@ill.fr>, September 2018.");
	tl2::log_info("--------------------------------------------------------------------------------");

	m_skx.SetT(-1000.);
	m_skx.SetB(25.);
	m_fp.SetT(29.);
	m_fp.SetB(0.35);
	m_heli.SetT(29.);
	m_heli.SetB(0.35);
	m_lf.SetT(29.);

	std::vector<ublas::vector<t_cplx>> fourier_skx;
	fourier_skx.reserve(_skxgs_allcomps.size()/3);

	for(std::size_t comp=0; comp<_skxgs_allcomps.size(); comp+=3)
		fourier_skx.push_back(tl2::make_vec<ublas::vector<t_cplx>>({_skxgs_allcomps[comp], _skxgs_allcomps[comp+1], _skxgs_allcomps[comp+2]}));

	m_skx.SetFourier(fourier_skx);
	m_skx.GenFullFourier();

	m_skx.SetCoords(m_vecB[0],m_vecB[1],m_vecB[2], m_vecPin[0],m_vecPin[1],m_vecPin[2]);
	m_skx.SetG(m_vecG[0], m_vecG[1], m_vecG[2]);
	m_skx.SetFilterZeroWeight(true);

	m_fp.SetCoords(m_vecB[0],m_vecB[1],m_vecB[2]);
	m_fp.SetG(m_vecG[0], m_vecG[1], m_vecG[2]);

	m_heli.SetCoords(m_vecB[0],m_vecB[1],m_vecB[2]);
	m_heli.SetG(m_vecG[0], m_vecG[1], m_vecG[2]);
	m_heli.SetFilterZeroWeight(true);

	m_lf.SetPinning(m_vecPin[0],m_vecPin[1],m_vecPin[2], m_vecB[0],m_vecB[1],m_vecB[2]);

	SqwBase::m_bOk = 1;
}

SqwMod::SqwMod(const std::string& strCfgFile) : SqwMod() { SqwBase::m_bOk = 1; }
SqwMod::~SqwMod() {}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// S(Q,E)

std::tuple<std::vector<t_real>, std::vector<t_real>>
	SqwMod::disp(t_real dh, t_real dk, t_real dl) const
{
	std::vector<t_real> vecE, vecW[4];
	if(m_iwhich_disp == 0)
		std::tie(vecE, vecW[0], vecW[1], vecW[2], vecW[3]) = m_skx.GetDisp(dh, dk, dl);
	else if(m_iwhich_disp == 1)
		std::tie(vecE, vecW[0], vecW[1], vecW[2], vecW[3]) = m_fp.GetDisp(dh, dk, dl);
	else if(m_iwhich_disp == 2)
		std::tie(vecE, vecW[0], vecW[1], vecW[2], vecW[3]) = m_heli.GetDisp(dh, dk, dl);

	return std::make_tuple(vecE, vecW[m_iPolChan]);
}


t_real SqwMod::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	t_vec vecq = tl2::make_vec<t_vec>({dh, dk, dl}) - m_vecG;

	// longitudinal fluctuations
	t_real dS_lf = m_lf.S_para(vecq, dE);
	if(m_ionlylf)
		return dS_lf;

	std::vector<t_real> vecE, vecW[4];
	// only calculate dispersion if global weight factor is not 0
	if(!tl2::float_equal(m_dS0, t_real(0)))
	{
		if(m_iwhich_disp == 0)
			std::tie(vecE, vecW[0], vecW[1], vecW[2], vecW[3]) = m_skx.GetDisp(dh, dk, dl, dE-m_dErange, dE+m_dErange);
		else if(m_iwhich_disp == 1)
			std::tie(vecE, vecW[0], vecW[1], vecW[2], vecW[3]) = m_fp.GetDisp(dh, dk, dl, dE-m_dErange, dE+m_dErange);
		else if(m_iwhich_disp == 2)
			std::tie(vecE, vecW[0], vecW[1], vecW[2], vecW[3]) = m_heli.GetDisp(dh, dk, dl, dE-m_dErange, dE+m_dErange);
	}

	// incoherent
	t_real dInc=0, dS_p=0, dS_m=0;
	if(!tl2::float_equal(m_dIncAmp, t_real(0)))
		dInc = tl2::gauss_model(dE, t_real(0), m_dIncSigma, m_dIncAmp, t_real(0));

	t_real dS = 0;
	for(std::size_t iE=0; iE<vecE.size(); ++iE)
	{
		if(!tl2::float_equal(vecW[m_iPolChan][iE], t_real(0)))
			dS += tl2::gauss_model(dE, vecE[iE], m_dSigma, vecW[m_iPolChan][iE], t_real(0));
	}

	return m_dS0*dS * tl2::bose_cutoff(dE, m_dT, m_dcut) + dInc;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// getters & setters

std::vector<SqwMod::t_var> SqwMod::GetVars() const
{
	std::vector<t_var> vecVars;

	vecVars.push_back(SqwBase::t_var{"proj_neutron", "real", tl2::var_to_str(m_iProjNeutron)});
	vecVars.push_back(SqwBase::t_var{"which_disp", "real", tl2::var_to_str(m_iwhich_disp)});
	vecVars.push_back(SqwBase::t_var{"E_range", "real", tl2::var_to_str(m_dErange)});
	vecVars.push_back(SqwBase::t_var{"bose_cutoff", "real", tl2::var_to_str(m_dcut)});
	vecVars.push_back(SqwBase::t_var{"only_mode", "int", tl2::var_to_str(m_iOnlyMode)});
	vecVars.push_back(SqwBase::t_var{"sigma", "real", tl2::var_to_str(m_dSigma)});
	vecVars.push_back(SqwBase::t_var{"inc_amp", "real", tl2::var_to_str(m_dIncAmp)});
	vecVars.push_back(SqwBase::t_var{"inc_sigma", "real", tl2::var_to_str(m_dIncSigma)});
	vecVars.push_back(SqwBase::t_var{"S0", "real", tl2::var_to_str(m_dS0)});
	vecVars.push_back(SqwBase::t_var{"T", "real", tl2::var_to_str(m_dT)});
	vecVars.push_back(SqwBase::t_var{"B", "real", tl2::var_to_str(m_dB)});
	vecVars.push_back(SqwBase::t_var{"pol_chan", "int", tl2::var_to_str(m_iPolChan)});
	vecVars.push_back(SqwBase::t_var{"G", "vector", vec_to_str(m_vecG)});
	vecVars.push_back(SqwBase::t_var{"B_dir", "vector", vec_to_str(m_vecB)});
	vecVars.push_back(SqwBase::t_var{"Pin_dir", "vector", vec_to_str(m_vecPin)});

	vecVars.push_back(SqwBase::t_var{"lf_only", "int", tl2::var_to_str(m_ionlylf)});
	vecVars.push_back(SqwBase::t_var{"lf_invcorrel", "real", tl2::var_to_str(m_lf.GetInvCorrel())});
	vecVars.push_back(SqwBase::t_var{"lf_A", "real", tl2::var_to_str(m_lf.GetA())});
	vecVars.push_back(SqwBase::t_var{"lf_gamma", "real", tl2::var_to_str(m_lf.GetGamma())});

	return vecVars;
}


void SqwMod::SetVars(const std::vector<SqwMod::t_var>& vecVars)
{
	if(!vecVars.size()) return;

	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "proj_neutron") m_iProjNeutron = tl2::str_to_var<int>(strVal);
		else if(strVar == "which_disp") m_iwhich_disp = tl2::str_to_var<int>(strVal);
		else if(strVar == "E_range") m_dErange = tl2::str_to_var<t_real>(strVal);
		else if(strVar == "bose_cutoff") m_dcut = tl2::str_to_var<t_real>(strVal);
		else if(strVar == "sigma") m_dSigma = tl2::str_to_var<t_real>(strVal);
		else if(strVar == "inc_amp") m_dIncAmp = tl2::str_to_var<decltype(m_dIncAmp)>(strVal);
		else if(strVar == "inc_sigma") m_dIncSigma = tl2::str_to_var<decltype(m_dIncSigma)>(strVal);
		else if(strVar == "S0") m_dS0 = tl2::str_to_var<decltype(m_dS0)>(strVal);
		else if(strVar == "only_mode")
		{
			m_iOnlyMode = tl2::str_to_var<int>(strVal);
			m_heli.SetOnlyMode(m_iOnlyMode);
		}
		else if(strVar == "T")
		{
			m_dT = tl2::str_to_var<decltype(m_dT)>(strVal);

			m_fp.SetT(m_dT);
			m_heli.SetT(m_dT);
			// fixed (and theo units) for skx!

			m_lf.SetT(m_dT);
		}
		else if(strVar == "B")
		{
			m_dB = tl2::str_to_var<decltype(m_dB)>(strVal);

			m_fp.SetB(m_dB);
			m_heli.SetB(m_dB);
			// fixed (and theo units) for skx!
		}
		else if(strVar == "pol_chan") m_iPolChan = tl2::str_to_var<decltype(m_iPolChan)>(strVal);
		else if(strVar == "G" || strVal == "proj_neutron")
		{
			m_vecG = str_to_vec<decltype(m_vecG)>(strVal);

			m_skx.SetProjNeutron(m_iProjNeutron!=0);
			m_skx.SetG(m_vecG[0], m_vecG[1], m_vecG[2]);

			m_fp.SetProjNeutron(m_iProjNeutron!=0);
			m_fp.SetG(m_vecG[0], m_vecG[1], m_vecG[2]);

			m_heli.SetProjNeutron(m_iProjNeutron!=0);
			m_heli.SetG(m_vecG[0], m_vecG[1], m_vecG[2]);
		}
		else if(strVar == "B_dir")
		{
			m_vecB = str_to_vec<decltype(m_vecB)>(strVal);

			m_skx.SetCoords(m_vecB[0],m_vecB[1],m_vecB[2], m_vecPin[0],m_vecPin[1],m_vecPin[2]);
			m_skx.SetG(m_vecG[0], m_vecG[1], m_vecG[2]);

			m_fp.SetCoords(m_vecB[0],m_vecB[1],m_vecB[2]);
			m_fp.SetG(m_vecG[0], m_vecG[1], m_vecG[2]);

			m_heli.SetCoords(m_vecB[0],m_vecB[1],m_vecB[2]);
			m_heli.SetG(m_vecG[0], m_vecG[1], m_vecG[2]);

			m_lf.SetPinning(m_vecPin[0],m_vecPin[1],m_vecPin[2],
				m_vecB[0],m_vecB[1],m_vecB[2]);
		}
		else if(strVar == "Pin_dir")
		{
			m_vecPin = str_to_vec<decltype(m_vecPin)>(strVal);

			m_skx.SetCoords(m_vecB[0],m_vecB[1],m_vecB[2], m_vecPin[0],m_vecPin[1],m_vecPin[2]);
			m_skx.SetG(m_vecG[0], m_vecG[1], m_vecG[2]);

			m_lf.SetPinning(m_vecPin[0],m_vecPin[1],m_vecPin[2],
				m_vecB[0],m_vecB[1],m_vecB[2]);
		}

		else if(strVar == "lf_only")
		{
			m_ionlylf = tl2::str_to_var<decltype(m_ionlylf)>(strVal);
		}
		else if(strVar == "lf_invcorrel")
		{
			t_real dInvCorrel = tl2::str_to_var<t_real>(strVal);
			m_lf.SetInvCorrel(dInvCorrel);
		}
		else if(strVar == "lf_A")
		{
			t_real dA = tl2::str_to_var<t_real>(strVal);
			m_lf.SetA(dA);
		}
		else if(strVar == "lf_gamma")
		{
			t_real dG = tl2::str_to_var<t_real>(strVal);
			m_lf.SetGamma(dG);
		}
	}
}

bool SqwMod::SetVarIfAvail(const std::string& strKey, const std::string& strNewVal)
{
	return SqwBase::SetVarIfAvail(strKey, strNewVal);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// copy

SqwBase* SqwMod::shallow_copy() const
{
	SqwMod *pMod = new SqwMod();

	pMod->m_iProjNeutron = this->m_iProjNeutron;
	pMod->m_iwhich_disp = this->m_iwhich_disp;
	pMod->m_dSigma = this->m_dSigma;
	pMod->m_dIncAmp = this->m_dIncAmp;
	pMod->m_dIncSigma = this->m_dIncSigma;
	pMod->m_dS0 = this->m_dS0;
	pMod->m_dT = this->m_dT;
	pMod->m_dB = this->m_dB;
	pMod->m_skx = this->m_skx;
	pMod->m_fp = this->m_fp;
	pMod->m_heli = this->m_heli;
	pMod->m_lf = this->m_lf;
	pMod->m_ionlylf = this->m_ionlylf;
	pMod->m_iPolChan = this->m_iPolChan;
	pMod->m_vecG = this->m_vecG;
	pMod->m_vecB = this->m_vecB;
	pMod->m_vecPin = this->m_vecPin;
	pMod->m_dErange = this->m_dErange;
	pMod->m_dcut = this->m_dcut;

	return pMod;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// takin interface

static const char* pcModIdent = "skxmod";
static const char* pcModName = "MnSi Magnon Dynamics";


#ifndef PLUGIN_APPLI


// exported symbols
std::tuple<std::string, std::string, std::string> sqw_info()
{
	return std::make_tuple(TAKIN_VER, pcModIdent, pcModName);
}

std::shared_ptr<SqwBase> sqw_construct(const std::string& strCfgFile)
{
	return std::make_shared<SqwMod>(strCfgFile);
}


// exports from so file
#ifndef __MINGW32__
	BOOST_DLL_ALIAS(sqw_info, takin_sqw_info);
	BOOST_DLL_ALIAS(sqw_construct, takin_sqw);
#else
	// hack because BOOST_DLL_ALIAS does not seem to work with Mingw

	extern "C" __declspec(dllexport)
	std::tuple<std::string, std::string, std::string> takin_sqw_info()
	{
		return sqw_info();
	}

	extern "C" __declspec(dllexport)
	std::shared_ptr<SqwBase> takin_sqw(const std::string& strCfgFile)
	{
		return sqw_construct(strCfgFile);
	}
#endif


#else


int main(int argc, char** argv)
{
	if(argc <= 2)
	{
		std::cout << "#\n# This is a Takin plugin module.\n#\n";
		std::cout << "module_ident: " << pcModIdent << "\n";
		std::cout << "module_name: " << pcModName << "\n";
		std::cout << "required_takin_version: " << TAKIN_VER << "\n";
		std::cout.flush();
		return 0;
	}

	const char* pcCfgFile = argv[1];
	const char* pcSharedMem = argv[2];
	SqwProc<SqwMod> proc(pcCfgFile, SqwProcStartMode::START_CHILD, pcSharedMem);
	return 0;
}


#endif
// ----------------------------------------------------------------------------
