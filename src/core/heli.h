/**
 * free energy and dynamics in the helimagnetic phase
 * @author tweber@ill.fr
 * @date mid-16, jul-18
 * @desc This file implements the theoretical helimagnon model by M. Garst and J. Waizner, references:
 *	- M. Garst, J. Waizner, and D. Grundler, J. Phys. D: Appl. Phys. 50, 293002 (2017), https://doi.org/10.1088/1361-6463/aa7573
 *	- J. Waizner, PhD thesis (2017), Universität zu Köln, https://kups.ub.uni-koeln.de/7937/
 *	- M. Kugler, G. Brandl, et al., Phys. Rev. Lett. 115, 097203 (2015), https://doi.org/10.1103/PhysRevLett.115.097203
 *	- Personal communications with M. Garst, 2014-2019.
 * @desc This file is based on:
 *	- The descriptions and Mathematica implementations of the different helimagnon model versions by M. Garst and J. Waizner, 2014-2018,
 *	- The 2015 and 2016 Python implementations by G. Brandl and M. Kugler of the first version of the helimagnon model.
 *	  This present version started as a C++ port of that Python implementation by G. Brandl and M. Kugler,
 *	  that was then adapted to new theoretical model revisions provided by M. Garst.
 *	- The 2016 Python implementations by M. Kugler and G. Brandl of the first version of the skyrmion model,
 *	  that also included improvements to the helimagnon code, which we ported to C++ in the present version.
 * @license GPLv2 (see 'LICENSE' file)
 */

#ifndef __HELI_H__
#define __HELI_H__

#include <tuple>
#include <vector>
#include <complex>

#include "helper.h"
#include "magsys.h"

#define USE_LAPACK
#include "tlibs2-extras/math17.h"


#ifndef DEF_HELI_ORDER
	#define DEF_HELI_ORDER 7
#endif


/**
 * helix
 */
template<class t_real=double, class t_cplx = std::complex<t_real>, int ORDER=4>
class Heli : public MagSystem<t_real, t_cplx, ORDER>, public MagDynamics<t_real, t_cplx>
{
public:
	static constexpr int ORDER_FOURIER = ORDER;
	using t_mat = ublas::matrix<t_real>;
	using t_vec = ublas::vector<t_real>;
	using t_quat = boost::math::quaternion<t_real>;
	using t_mat_cplx = ublas::matrix<t_cplx>;
	using t_vec_cplx = ublas::vector<t_cplx>;


public:
	Heli();
	virtual ~Heli() = default;

	virtual std::shared_ptr<Heli<t_real, t_cplx, ORDER>> copy() const
	{
		return std::make_shared<Heli<t_real, t_cplx, ORDER>>(*this);
	}

	virtual std::shared_ptr<MagSystem<t_real, t_cplx, ORDER_FOURIER>> copyCastSys() const override
	{
		return copy();
	}

	virtual std::shared_ptr<MagDynamics<t_real, t_cplx>> copyCastDyn() const override
	{
		return copy();
	}

	virtual t_real F() override;

	virtual void SetFourier(const std::vector<ublas::vector<t_cplx>> &fourier) override;
	virtual const std::vector<ublas::vector<t_cplx>> &GetFourier() const override { return m_fourier; }

	// careful: free energy still uses theory units for T and B, dynamics calculation already uses real units!
	virtual void SetB(t_real B) override { m_B = B; }
	virtual void SetT(t_real T) override { m_T = T;  m_Bc2 = get_bc2(m_T, false); }
	t_real GetBC2() const { return m_Bc2; }

	using MagSystem<t_real, t_cplx, ORDER_FOURIER>::minimise;


	void SetG(t_real h, t_real k, t_real l);
	void SetCoords(t_real Bx, t_real By, t_real Bz);
	void SetFilterZeroWeight(bool b) { m_filterzeroweight = b; }
	void SetWeightEps(t_real eps) { m_weighteps = eps; }
	void SetOnlyMode(int iMode) { m_onlymode = iMode; }
	void SetProjNeutron(bool b) { m_bProjNeutron = b; }

	std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
	GetSpecWeights(t_real qh, t_real qk, t_real ql, t_real minE=-1., t_real maxE=-2.) const;

	virtual std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
		GetDisp(t_real h, t_real k, t_real l, t_real minE=-1., t_real maxE=-2.) const override;


private:
	t_real m_pitch = 1;
	t_real m_B = 0, m_T = -100;

	std::vector<ublas::vector<t_cplx>> m_fourier{};


	int m_onlymode = -1;
	bool m_filterzeroweight = false;
	t_real m_eveps = 1e-6;
	t_real m_weighteps = 1e-6;

	t_real m_Bc2 = 0;
	t_vec m_Grlu{};
	t_quat m_rotCoord{};

	std::vector<t_mat_cplx> m_polMat =
	{{
		get_polmat<t_mat_cplx>(2),	// SF1
		get_polmat<t_mat_cplx>(1),	// SF2
		get_polmat<t_mat_cplx>(3)	// NSF
	}};
	t_mat_cplx m_projNeutron = tl2::unit_m<t_mat_cplx>(3);
	bool m_bProjNeutron = true;


private:
	std::vector<int> m_idx2[3], m_idx3[4];
};


#endif
