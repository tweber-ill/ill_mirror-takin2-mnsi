/**
 * Takin module for precalculated grid data
 * @author Tobias Weber <tweber@ill.fr>
 * @date nov-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "tools/monteconvo/sqwbase.h"
#include "libs/version.h"

#include "tlibs2/libs/math17.h"
#include "tlibs2/libs/str.h"
#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/file.h"

#include <boost/dll/alias.hpp>
#include <QtCore/QFile>


using t_real = t_real_reso;
using t_vec = ublas::vector<t_real>;
using t_mat = ublas::matrix<t_real>;
using t_quat = boost::math::quaternion<t_real>;
using t_cplx = std::complex<t_real>;


class SqwMod : public SqwBase
{
	public:
		using SqwBase::t_var;

	protected:
		t_real m_dSigma = t_real(0.05);
		t_real m_dS0 = t_real(1.);
		t_real m_dIncAmp = t_real(0.);
		t_real m_dIncSigma = t_real(0.05);
		t_real m_dT = 29.;
		t_real m_dRotZ = 0.;
		t_real m_dRotX = 0.;
		t_real m_dcut = 0.02;
		int m_iPolChan = 0;
		bool m_bFlipB_or_q = false;
		t_real m_dEFactor = 0.24;
		t_real m_dEFactorGrid = 0.24;

		t_vec m_vecG = tl2::make_vec<t_vec>({1,1,0});

		t_vec m_vecRotFrom = tl2::make_vec<t_vec>({1,0,0});
		t_vec m_vecRotTo = tl2::make_vec<t_vec>({1,0,0});

		// grid descriptors
		std::string m_strIndexFile, m_strDataFile;

		t_real m_hmin=0., m_hmax=0., m_hstep=0.;
		t_real m_kmin=0., m_kmax=0., m_kstep=0.;
		t_real m_lmin=0., m_lmax=0., m_lstep=0.;


	public:
		SqwMod();
		SqwMod(const std::string& strCfgFile);

		virtual std::tuple<std::vector<t_real>, std::vector<t_real>>
			disp(t_real dh, t_real dk, t_real dl) const override;
		virtual t_real operator()(t_real dh, t_real dk, t_real dl, t_real dE) const override;

		virtual std::vector<t_var> GetVars() const override;
		virtual void SetVars(const std::vector<t_var>&) override;
		virtual bool SetVarIfAvail(const std::string& strKey, const std::string& strNewVal) override;

		virtual SqwBase* shallow_copy() const override;
};



// ----------------------------------------------------------------------------
// constructors

SqwMod::SqwMod()
{
	SqwBase::m_bOk = 1;
}


SqwMod::SqwMod(const std::string& strCfgFile) : SqwMod()
{
	tl2::log_info("--------------------------------------------------------------------------------");
	tl2::log_info("This is the Takin skx grid module,");
	tl2::log_info("\tby T. Weber <tweber@ill.fr>, November 2018.");
	tl2::log_info("--------------------------------------------------------------------------------");

	tl2::log_info("Grid description file: \"", strCfgFile, "\".");

	tl2::Prop<std::string> prop;
	if(prop.Load(strCfgFile.c_str(), tl2::PropType::INFO))
	{
		m_strIndexFile = prop.Query<std::string>("files/index");
		m_strDataFile = prop.Query<std::string>("files/data");

		m_hmin = prop.QueryAndParse<t_real>("dims/hmin");
		m_hmax = prop.QueryAndParse<t_real>("dims/hmax");
		m_hstep = prop.QueryAndParse<t_real>("dims/hstep");
		m_kmin = prop.QueryAndParse<t_real>("dims/kmin");
		m_kmax = prop.QueryAndParse<t_real>("dims/kmax");
		m_kstep = prop.QueryAndParse<t_real>("dims/kstep");
		m_lmin = prop.QueryAndParse<t_real>("dims/lmin");
		m_lmax = prop.QueryAndParse<t_real>("dims/lmax");
		m_lstep = prop.QueryAndParse<t_real>("dims/lstep");

		SqwBase::m_bOk = 1;

		tl2::log_info("Index file: ", m_strIndexFile, ", data file: ", m_strDataFile, ".");
		tl2::log_info("h range: ", m_hmin, " .. ", m_hmax, ", step: ", m_hstep);
		tl2::log_info("k range: ", m_kmin, " .. ", m_kmax, ", step: ", m_kstep);
		tl2::log_info("l range: ", m_lmin, " .. ", m_lmax, ", step: ", m_lstep);
	}
	else
	{
		tl2::log_err("Grid description file \"",  strCfgFile, "\" could not be loaded.");
		SqwBase::m_bOk = 0;
	}
}



// ----------------------------------------------------------------------------
// dispersion, spectral weight and structure factor

std::tuple<std::vector<t_real>, std::vector<t_real>>
	SqwMod::disp(t_real dh, t_real dk, t_real dl) const
{
	if(!tl2::vec_equal(m_vecRotFrom, m_vecRotTo, 1e-6))
	{
		t_quat quat = tl2::rotation_quat(m_vecRotFrom, m_vecRotTo);
		t_vec vecPos = tl2::make_vec<t_vec>({dh, dk, dl}) - m_vecG;
		vecPos = tl2::quat_vec_prod(quat, vecPos) + m_vecG;
		dh = vecPos[0];
		dk = vecPos[1];
		dl = vecPos[2];
	}

	if(!tl2::float_equal(m_dRotZ, 0., 1e-6))
	{
		t_mat rot = tl2::rotation_matrix_3d_z<t_mat>(m_dRotZ/180.*tl2::pi<t_real>);
		t_vec vecPos = tl2::make_vec<t_vec>({dh, dk, dl}) - m_vecG;
		vecPos = tl2::prod_mv(rot, vecPos) + m_vecG;

		dh = vecPos[0];
		dk = vecPos[1];
		dl = vecPos[2];
	}

	if(!tl2::float_equal(m_dRotX, 0., 1e-6))
	{
		t_mat rot = tl2::rotation_matrix_3d_x<t_mat>(m_dRotZ/180.*tl2::pi<t_real>);
		t_vec vecPos = tl2::make_vec<t_vec>({dh, dk, dl}) - m_vecG;
		vecPos = tl2::prod_mv(rot, vecPos) + m_vecG;

		dh = vecPos[0];
		dk = vecPos[1];
		dl = vecPos[2];
	}


	/**
	 * calculate file index based on coordinates
	 */
	auto hkl_to_idx = [this](t_real h, t_real k, t_real l) -> std::size_t
	{
		// clamp values to boundaries
		if(h < m_hmin) h = m_hmin;
		if(k < m_kmin) k = m_kmin;
		if(l < m_lmin) l = m_lmin;
		if(h >= m_hmax) h = m_hmax - m_hstep;
		if(k >= m_kmax) k = m_kmax - m_kstep;
		if(l >= m_lmax) l = m_lmax - m_lstep;

		// max dimensions
		std::size_t iHSize = std::size_t(((m_hmax-m_hmin) / m_hstep));
		std::size_t iKSize = std::size_t(((m_kmax-m_kmin) / m_kstep));
		std::size_t iLSize = std::size_t(((m_lmax-m_lmin) / m_lstep));

		// position indices
		std::size_t iH = std::size_t(std::round(((h - m_hmin) / m_hstep)));
		std::size_t iK = std::size_t(std::round(((k - m_kmin) / m_kstep)));
		std::size_t iL = std::size_t(std::round(((l - m_lmin) / m_lstep)));

		// clamp again
		if(iH >= iHSize) iH = iHSize-1;
		if(iK >= iKSize) iK = iKSize-1;
		if(iL >= iLSize) iL = iLSize-1;

		return iH*iKSize*iLSize + iK*iLSize + iL;
	};


	// ------------------------------------------------------------------------
	// the index file holds the offsets into the data file
	t_real dqh = dh-m_vecG[0], dqk = dk-m_vecG[1], dql = dl-m_vecG[2];
	if(m_bFlipB_or_q)
	{
		dqh = -dqh;
		dqk = -dqk;
		dql = -dql;
	}
	std::size_t idx_file_offs = hkl_to_idx(dqh, dqk, dql);

	QFile fileIdx(m_strIndexFile.c_str());
	if(!fileIdx.exists())
	{
		tl2::log_err("Index file \"", m_strIndexFile, "\" does not exist.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	if(!fileIdx.open(QIODevice::ReadOnly))
	{
		tl2::log_err("Index file \"", m_strIndexFile, "\" cannot be opened.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	const void *pMemIdx = fileIdx.map(idx_file_offs*sizeof(std::size_t), sizeof(std::size_t));
	if(!pMemIdx)
	{
		tl2::log_err("Index file \"", m_strIndexFile, "\" cannot be mapped.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	std::size_t dat_file_offs = *((std::size_t*)pMemIdx);
	fileIdx.unmap((unsigned char*)pMemIdx);
	fileIdx.close();
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// the data file holds the energies and spectral weights of the dispersion branches

	QFile fileDat(m_strDataFile.c_str());
	if(!fileDat.exists())
	{
		tl2::log_err("Data file \"", m_strDataFile, "\" does not exist.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	if(!fileDat.open(QIODevice::ReadOnly))
	{
		tl2::log_err("Data file \"", m_strDataFile, "\" cannot be opened.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	const void *pMemDat = fileDat.map(dat_file_offs, sizeof(std::size_t));
	if(!pMemDat)
	{
		tl2::log_err("Data file \"", m_strDataFile, "\" cannot be mapped (1).");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	// number of dispersion branches and weights
	unsigned int iNumBranches = *((unsigned int*)pMemDat);
	fileDat.unmap((unsigned char*)pMemDat);


	// map actual (E, w) data
	const t_real *pBranches = (t_real*)fileDat.map(dat_file_offs+sizeof(iNumBranches), iNumBranches*sizeof(t_real)*4);
	if(!pMemDat)
	{
		tl2::log_err("Data file \"", m_strDataFile, "\" cannot be mapped (2).");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}


	std::vector<t_real> vecE, vecw, vecwSF1, vecwSF2, vecwNSF;
	for(unsigned int iBranch=0; iBranch<iNumBranches; ++iBranch)
	{
		vecE.push_back(pBranches[iBranch*4 + 0] * m_dEFactor / m_dEFactorGrid);	// energy
		vecwSF1.push_back(pBranches[iBranch*4 + 1]);	// weight SF1
		vecwSF2.push_back(pBranches[iBranch*4 + 2]);	// weight SF2
		vecwNSF.push_back(pBranches[iBranch*4 + 3]);	// weight NSF
		vecw.push_back(pBranches[iBranch*4 + 1] + pBranches[iBranch*4 + 2] + pBranches[iBranch*4 + 3]);	// full weight
	}

	fileDat.unmap((unsigned char*)pBranches);
	fileDat.close();
	// ------------------------------------------------------------------------


	if(m_iPolChan == 1) return std::make_tuple(vecE, vecwSF1);
	else if(m_iPolChan == 2) return std::make_tuple(vecE, vecwSF2);
	else if(m_iPolChan == 3) return std::make_tuple(vecE, vecwNSF);
	return std::make_tuple(vecE, vecw);
}


/**
 * S(q,E)
 */
t_real SqwMod::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	t_real dInc=0, dS_p=0, dS_m=0;
	if(!tl2::float_equal(m_dIncAmp, t_real(0)))
		dInc = tl2::gauss_model(dE, t_real(0), m_dIncSigma, m_dIncAmp, t_real(0));

	t_real dS = 0;
	if(!tl2::float_equal(m_dS0, t_real(0)))
	{
		std::vector<t_real> vecE, vecW;
		std::tie(vecE, vecW) = disp(dh, dk, dl);

		for(std::size_t iE=0; iE<vecE.size(); ++iE)
			dS += tl2::gauss_model(dE, vecE[iE], m_dSigma, vecW[iE], t_real(0));
	}

	return m_dS0*dS * tl2::bose_cutoff(dE, m_dT, m_dcut) + dInc;
}



// ----------------------------------------------------------------------------
// get & set variables

std::vector<SqwMod::t_var> SqwMod::GetVars() const
{
	std::vector<t_var> vecVars;

	vecVars.push_back(SqwBase::t_var{"bose_cutoff", "real", tl2::var_to_str(m_dcut)});
	vecVars.push_back(SqwBase::t_var{"sigma", "real", tl2::var_to_str(m_dSigma)});
	vecVars.push_back(SqwBase::t_var{"inc_amp", "real", tl2::var_to_str(m_dIncAmp)});
	vecVars.push_back(SqwBase::t_var{"inc_sigma", "real", tl2::var_to_str(m_dIncSigma)});
	vecVars.push_back(SqwBase::t_var{"S0", "real", tl2::var_to_str(m_dS0)});
	vecVars.push_back(SqwBase::t_var{"T", "real", tl2::var_to_str(m_dT)});
	vecVars.push_back(SqwBase::t_var{"rot_z", "real", tl2::var_to_str(m_dRotZ)});
	vecVars.push_back(SqwBase::t_var{"rot_x", "real", tl2::var_to_str(m_dRotX)});
	vecVars.push_back(SqwBase::t_var{"pol_chan", "int", tl2::var_to_str(m_iPolChan)});
	vecVars.push_back(SqwBase::t_var{"G", "vector", vec_to_str(m_vecG)});
	vecVars.push_back(SqwBase::t_var{"rot_from", "vector", vec_to_str(m_vecRotFrom)});
	vecVars.push_back(SqwBase::t_var{"rot_to", "vector", vec_to_str(m_vecRotTo)});
	vecVars.push_back(SqwBase::t_var{"skx_E_factor", "real", tl2::var_to_str(m_dEFactor)});
	vecVars.push_back(SqwBase::t_var{"flip_B_or_q", "int", tl2::var_to_str(int(m_bFlipB_or_q))});

	return vecVars;
}


void SqwMod::SetVars(const std::vector<SqwMod::t_var>& vecVars)
{
	if(!vecVars.size()) return;

	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "bose_cutoff") m_dcut = tl2::str_to_var<t_real>(strVal);
		else if(strVar == "sigma") m_dSigma = tl2::str_to_var<t_real>(strVal);
		else if(strVar == "inc_amp") m_dIncAmp = tl2::str_to_var<decltype(m_dIncAmp)>(strVal);
		else if(strVar == "inc_sigma") m_dIncSigma = tl2::str_to_var<decltype(m_dIncSigma)>(strVal);
		else if(strVar == "S0") m_dS0 = tl2::str_to_var<decltype(m_dS0)>(strVal);
		else if(strVar == "T") m_dT = tl2::str_to_var<decltype(m_dT)>(strVal);
		else if(strVar == "rot_z") m_dRotZ = tl2::str_to_var<decltype(m_dRotZ)>(strVal);
		else if(strVar == "rot_x") m_dRotZ = tl2::str_to_var<decltype(m_dRotX)>(strVal);
		else if(strVar == "pol_chan") m_iPolChan = tl2::str_to_var<decltype(m_iPolChan)>(strVal);
		else if(strVar == "G") m_vecG = str_to_vec<decltype(m_vecG)>(strVal);
		else if(strVar == "rot_from") m_vecRotFrom = str_to_vec<decltype(m_vecRotFrom)>(strVal);
		else if(strVar == "rot_to") m_vecRotTo = str_to_vec<decltype(m_vecRotTo)>(strVal);
		else if(strVar == "skx_E_factor") m_dEFactor = tl2::str_to_var<decltype(m_dT)>(strVal);
		else if(strVar == "flip_B_or_q") m_bFlipB_or_q = (tl2::str_to_var<int>(strVal) != 0);
	}
}


bool SqwMod::SetVarIfAvail(const std::string& strKey, const std::string& strNewVal)
{
	return SqwBase::SetVarIfAvail(strKey, strNewVal);
}



// ----------------------------------------------------------------------------
// copy

SqwBase* SqwMod::shallow_copy() const
{
	SqwMod *pMod = new SqwMod();

	pMod->m_dSigma = this->m_dSigma;
	pMod->m_dIncAmp = this->m_dIncAmp;
	pMod->m_dIncSigma = this->m_dIncSigma;
	pMod->m_dS0 = this->m_dS0;
	pMod->m_dT = this->m_dT;
	pMod->m_dRotZ = this->m_dRotZ;
	pMod->m_dRotX = this->m_dRotX;
	pMod->m_iPolChan = this->m_iPolChan;
	pMod->m_vecG = this->m_vecG;
	pMod->m_vecRotFrom = this->m_vecRotFrom;
	pMod->m_vecRotTo = this->m_vecRotTo;
	pMod->m_dcut = this->m_dcut;
	pMod->m_dEFactor = this->m_dEFactor;
	pMod->m_dEFactorGrid = this->m_dEFactorGrid;
	pMod->m_bFlipB_or_q = this->m_bFlipB_or_q;

	pMod->m_strIndexFile = this->m_strIndexFile;
	pMod->m_strDataFile = this->m_strDataFile;

	pMod->m_hmin = this->m_hmin;
	pMod->m_hmax = this->m_hmax;
	pMod->m_hstep = this->m_hstep;
	pMod->m_kmin = this->m_kmin;
	pMod->m_kmax = this->m_kmax;
	pMod->m_kstep = this->m_kstep;
	pMod->m_lmin = this->m_lmin;
	pMod->m_lmax = this->m_lmax;
	pMod->m_lstep = this->m_lstep;

	return pMod;
}



#ifdef DO_TEST
// ----------------------------------------------------------------------------
// test querying of data

int main(int argc, char **argv)
{
	if(argc <= 1)
	{
		std::cerr << "Please specify a config file." << std::endl;
		return -1;
	}


	SqwMod mod(argv[1]);
	mod.SetVarIfAvail("pol_chan", "0");

	while(1)
	{
		t_real h, k, l;
		std::cout << "Enter hkl: ";
		std::cin >> h >> k >> l;

		std::vector<t_real> vecE, vecW;
		std::tie(vecE, vecW) = mod.disp(h, k, l);

		for(std::size_t i=0; i<vecE.size(); ++i)
			std::cout << "(" << i+1 << ") " << "E = " << vecE[i] << ", weight = " << vecW[i] << "\n";
		std::cout << std::endl;
	}
}

// ----------------------------------------------------------------------------
#else
// ----------------------------------------------------------------------------
// SO interface

static const char* pcModIdent = "skxmod_grid";
static const char* pcModName = "MnSi Magnon Dynamics (Grid)";


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
// ----------------------------------------------------------------------------


#endif
