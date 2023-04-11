/**
 * free energy and dynamics in the skyrmion phase
 * @author tweber@ill.fr
 * @date jul-18
 * @desc This file implements the theoretical skyrmion model by M. Garst and J. Waizner, references:
 *	- M. Garst, J. Waizner, and D. Grundler, J. Phys. D: Appl. Phys. 50, 293002 (2017), https://doi.org/10.1088/1361-6463/aa7573
 *	- J. Waizner, PhD thesis (2017), Universität zu Köln, https://kups.ub.uni-koeln.de/7937/
 *	- Personal communications with M. Garst, 2017-2020.
 * @desc This file is based on:
 *	- The descriptions and Mathematica implementations of the different skyrmion model versions by M. Garst and J. Waizner, 2016-2020.
 *	- The 2016 optimised Python implementations by M. Kugler and G. Brandl of the first version of the skyrmion model.
 *	  This present version started as a C++ port of M. Kugler's and G. Brandl's Python implementation,
 *	  that was then adapted to new theoretical model revisions provided by M. Garst.
 * @license GPLv2 (see 'LICENSE' file)
 */

#ifndef __SKX_H__
#define __SKX_H__

#include "magsys.h"
#define USE_LAPACK
#include "tlibs2/libs/math17.h"


#ifndef DEF_SKX_ORDER
	#define DEF_SKX_ORDER 7
#endif

#ifndef SKX_USE_HOC
	#define SKX_USE_HOC 0
#endif


/**
 * skyrmion dynamics
 */
template<class t_real=double, class t_cplx = std::complex<t_real>, int ORDER=4>
class Skx : public MagSystem<t_real, t_cplx, (ORDER+1)*ORDER/2>, public MagDynamics<t_real, t_cplx>
{
public:
	static constexpr int SIZE = 2*ORDER + 1;                // size of [-ORDER, ORDER] range
	static constexpr int ORDER_FOURIER = (ORDER+1)*ORDER/2; // number of peaks in 60 degree segment

	using t_mat = ublas::matrix<t_real>;
	using t_vec = ublas::vector<t_real>;
	using t_quat = boost::math::quaternion<t_real>;
	using t_mat_cplx = ublas::matrix<t_cplx>;
	using t_vec_cplx = ublas::vector<t_cplx>;
	using t_vec_int = ublas::vector<int>;


public:
	Skx();
	virtual ~Skx() = default;

	virtual t_real F() override;

	virtual void SetFourier(const std::vector<t_vec_cplx> &fourier, bool symm=true) override;
	virtual const std::vector<t_vec_cplx>& GetFourier() const override { return m_fourier; }

	virtual void SetB(t_real B, bool exp=true) override
	{
		if(!exp)
			m_B = B;
	}

	virtual void SetT(t_real T, bool exp=true) override
	{
		if(exp)
		{
			m_Bc2_exp = get_bc2(T, !exp, !SKX_USE_HOC);
		}
		else
		{
			m_T = T;
			m_Bc2 = get_bc2(m_T, !exp, !SKX_USE_HOC);
		}
	}

	virtual t_real GetBC2(bool exp=true) const override
	{
		return exp ? m_Bc2_exp : m_Bc2;
	}

	using MagSystem<t_real, t_cplx, ORDER_FOURIER>::minimise;


	std::shared_ptr<Skx<t_real, t_cplx, ORDER>> copy() const
	{
		return std::make_shared<Skx<t_real, t_cplx, ORDER>>(*this);
	}

	virtual std::shared_ptr<MagSystem<t_real, t_cplx, ORDER_FOURIER>> copyCastSys() const override
	{
		return copy();
	}

	virtual std::shared_ptr<MagDynamics<t_real, t_cplx>> copyCastDyn() const override
	{
		return copy();
	}


	virtual void SetCoords(t_real Bx, t_real By, t_real Bz, t_real Pinx, t_real Piny, t_real Pinz) override;

	virtual void SetG(t_real h, t_real k, t_real l) override
	{
		m_Grlu = tl2::make_vec<t_vec>({ h, k, l });
	}

	virtual void SetProjNeutron(bool b) override { m_bProjNeutron = b; }
	virtual void SetFilterZeroWeight(bool b) override { m_filterzeroweight = b; }

	virtual std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
		GetDisp(t_real h, t_real k, t_real l, t_real minE=-1., t_real maxE=-2.) const override;

	void SetWeightEps(t_real eps) { m_weighteps = eps; }

	const std::vector<t_vec>& GetPeaks60() { return m_peaks60_rlu; }


protected:
	std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
		GetSpecWeights(int Ghmag, int Gkmag, t_real qh, t_real qk, t_real ql,
			const t_mat_cplx& projNeutron,
			t_real minE=-1., t_real maxE=-2.) const;


private:
	bool m_filterzeroweight = false;
	t_real m_eps = 1e-5;
	t_real m_eveps = 1e-6;
	t_real m_weighteps = 1e-6;
	t_real m_evlimit = 0.9995;

	// B matrix with 120 deg between (100) and (010), Q_lab = B*Q_rlu
	t_mat m_Bmat = tl2::make_mat<t_mat>({
		{ 1, std::cos(tl2::d2r<t_real>(120)) },
		{ 0, std::sin(tl2::d2r<t_real>(120)) }
	});
	t_mat m_Binv{};

	t_real m_B = 0;
	t_real m_T = -1000;
	t_real m_Bc2 = 0;
	t_real m_Bc2_exp = 0.322;

	std::vector<t_vec_cplx> m_fourier{};
	std::vector<std::pair<int, int>> m_idx2[3], m_idx3[4];

	// all satellite peaks and one sixth of them, respectively
	std::vector<t_vec> m_allpeaks_rlu{}, m_peaks60_rlu{};

	ublas::matrix<t_vec_cplx> m_M{};
	t_vec m_Grlu{};
	t_quat m_rotCoord{};

	std::vector<t_mat_cplx> m_polMat =
	{{
		get_chiralpol<t_mat_cplx>(1),	// SF1
		get_chiralpol<t_mat_cplx>(2),	// SF2
		get_chiralpol<t_mat_cplx>(3)	// NSF
	}};

	bool m_bProjNeutron = true;
};


#endif
