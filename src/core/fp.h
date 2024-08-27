/**
 * dynamics in the field-polarised phase
 * @author Tobias Weber <tweber@ill.fr>
 * @date dec-16, sep-18
 * @desc This file implements the theoretical magnon model by M. Garst and J. Waizner, see:
 *	- M. Garst, J. Waizner, and D. Grundler, J. Phys. D: Appl. Phys. 50, 293002 (2017), https://doi.org/10.1088/1361-6463/aa7573
 *	- J. Waizner, PhD thesis (2017), Universität zu Köln, https://kups.ub.uni-koeln.de/7937/
 *	- personal communications with M. Garst, 2016-2019.
 * @desc This file is based on:
 *	- the descriptions and Mathematica implementations of the field-polarised magnon model by M. Garst, 2016.
 * @license GPLv2 (see 'LICENSE' file)
 */

#ifndef __FP_H__
#define __FP_H__

#include <tuple>
#include <vector>
#include <complex>

#include "helper.h"
#include "magsys.h"

#define USE_LAPACK
#include "tlibs2/libs/math17.h"


/**
 * field-polarised dynamics
 */
template<class t_real=double, class t_cplx = std::complex<t_real>>
class FP : public MagDynamics<t_real, t_cplx>
{
public:
	using t_mat = ublas::matrix<t_real>;
	using t_vec = ublas::vector<t_real>;
	using t_quat = boost::math::quaternion<t_real>;
	using t_mat_cplx = ublas::matrix<t_cplx>;
	using t_vec_cplx = ublas::vector<t_cplx>;


public:
	FP();
	virtual ~FP() = default;

	std::shared_ptr<FP<t_real, t_cplx>> copy() const
	{
		return std::make_shared<FP<t_real, t_cplx>>(*this);
	}

	virtual std::shared_ptr<MagDynamics<t_real, t_cplx>> copyCastDyn() const override
	{
		return copy();
	}

	virtual void SetB(t_real B, bool exp = true) override
	{
		if(exp)
			m_B = B;
	}

	virtual void SetT(t_real T, bool exp = true) override
	{
		if(exp)
		{
			m_T = T;
			m_Bc2 = get_bc2(m_T, false);
		}
	}

	virtual t_real GetBC2(bool exp=true) const override
	{
		return exp ? m_Bc2 : 0.;
	}

	virtual void SetG(t_real h, t_real k, t_real l) override
	{
		m_Grlu = tl2::make_vec<t_vec>({ h, k, l });
	}

	virtual void SetCoords(t_real Bx, t_real By, t_real Bz, t_real Px=0., t_real Py=0., t_real Pz=0.) override;
	virtual void SetProjNeutron(bool b) override { m_bProjNeutron = b; }
	virtual void SetFilterZeroWeight(bool) override {}

	virtual std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
		GetDisp(t_real h, t_real k, t_real l, t_real minE=-1., t_real maxE=-2.) const override;


private:
	// constants
	t_real m_Bc2 = 0.55;
	t_real m_B = 0.7;
	t_real m_T = 20;

	t_vec m_Grlu{};
	t_quat m_rotCoord{};

	bool m_bProjNeutron = true;

	static const t_vec_cplx m_x, m_y;
	static const t_vec m_z;
	static const t_mat_cplx m_sigma2, m_ident2;
};


#endif
