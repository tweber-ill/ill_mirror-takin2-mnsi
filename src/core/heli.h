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
 *	- The 2015 and 2016 optimised Python implementations by G. Brandl and M. Kugler of the first version of the helimagnon model.
 *	  This present version started as a C++ port of G. Brandl's and M. Kugler's Python implementation,
 *	  that was then adapted to new theoretical model revisions provided by M. Garst.
 *	- The 2016 optimised Python implementations by M. Kugler and G. Brandl of the first version of the skyrmion model,
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
#include "tlibs2/libs/math17.h"

#ifndef DEF_HELI_ORDER
	#define DEF_HELI_ORDER 7
#endif

#ifndef HELI_USE_HOC
	#define HELI_USE_HOC 0
#endif


/**
 * helical dynamics
 */
template<class t_real=double, class t_cplx = std::complex<t_real>, int ORDER=4>
class Heli : public MagSystem<t_real, t_cplx, ORDER>, public MagDynamics<t_real, t_cplx>
{
public:
	static constexpr int SIZE = 2*ORDER+1; // size of [-ORDER, ORDER] range
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

	virtual void SetFourier(const std::vector<t_vec_cplx> &fourier, bool symm=true) override;
	virtual const std::vector<t_vec_cplx>& GetFourier() const override { return m_fourier; }

	virtual void SetB(t_real B, bool exp=true) override
	{
		if(exp)
			m_B = B;
		else
			m_B_theo = B;
	}

	virtual void SetT(t_real T, bool exp=true) override
	{
		if(exp)
		{
			m_T = T;
			m_Bc2 = get_bc2(m_T, !exp, !HELI_USE_HOC);
		}
		else
		{
			m_T_theo = T;
			m_Bc2_theo = get_bc2(m_T_theo, !exp, !HELI_USE_HOC);
		}
	}

	virtual t_real GetBC2(bool exp=true) const override
	{
		return exp ? m_Bc2 : m_Bc2_theo;
	}

	using MagSystem<t_real, t_cplx, ORDER_FOURIER>::minimise;


	virtual void SetG(t_real h, t_real k, t_real l) override
	{
		m_Grlu = tl2::make_vec<t_vec>({ h, k, l });
	}

	virtual void SetCoords(t_real Bx, t_real By, t_real Bz, t_real Px=0., t_real Py=0., t_real Pz=0.) override;
	virtual void SetFilterZeroWeight(bool b) override { m_filterzeroweight = b; }
	virtual void SetProjNeutron(bool b) override { m_projNeutron = b; }

	bool GetExplicitCalc() const { return m_explicitcalc; }
	void SetExplicitCalc(bool b) { m_explicitcalc = b; }

	virtual std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
		GetDisp(t_real h, t_real k, t_real l, t_real minE=-1., t_real maxE=-2.) const override;

	void SetWeightEps(t_real eps) { m_weighteps = eps; }
	void SetOnlyMode(int iMode) { m_onlymode = iMode; }

	std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
		GetSpecWeights(t_real qh, t_real qk, t_real ql,
			const t_mat_cplx& projNeutron,
			t_real minE = -1., t_real maxE = -2.) const;


private:
	t_real m_eps = 1e-5;
	t_real m_B = 0.2, m_B_theo = 25.;
	t_real m_T = 27.5, m_T_theo = -1000.;
	t_real m_Bc2 = 0.3, m_Bc2_theo = 45.;

	std::vector<t_vec_cplx> m_fourier{};
	std::vector<int> m_idx2[3], m_idx3[4];

	int m_onlymode = -1;
	bool m_filterzeroweight = false;
	t_real m_eveps = 1e-6;
	t_real m_weighteps = 1e-6;

	t_vec m_Grlu{};
	t_quat m_rotCoord{};

	bool m_projNeutron = true;
	bool m_explicitcalc = true;
};


#endif
