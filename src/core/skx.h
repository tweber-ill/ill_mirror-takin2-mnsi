/**
 * Free energy in the skyrmion phase
 * @author tweber@ill.fr
 * @date jul-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#ifndef __SKX_H__
#define __SKX_H__

#include "magsys.h"
#define USE_LAPACK
//#include "tlibs2/libs/math17.h"
#include "tlibs2-extras/math17.h"


#ifndef DEF_SKX_ORDER
	#define DEF_SKX_ORDER 7
#endif


/**
 * skyrmion
 */
template<class t_real=double, class t_cplx = std::complex<t_real>, int ORDER=4>
class Skx : public MagSystem<t_real, t_cplx, (ORDER+1)*ORDER/2>, public MagDynamics<t_real, t_cplx>
{
public:
	static constexpr int ORDER_FOURIER = (ORDER+1)*ORDER/2;

	using t_mat = ublas::matrix<t_real>;
	using t_vec = ublas::vector<t_real>;
	using t_quat = boost::math::quaternion<t_real>;
	using t_mat_cplx = ublas::matrix<t_cplx>;
	using t_vec_cplx = ublas::vector<t_cplx>;


public:
	Skx();

	virtual t_real F() override;

	virtual void SetFourier(const std::vector<t_vec_cplx> &fourier) override;
	virtual const std::vector<t_vec_cplx> &GetFourier() const override { return m_fourier; }

	virtual void SetB(t_real B) override { m_B = B; }
	virtual void SetT(t_real T) override { m_T = T; m_Bc2 = get_bc2(m_T); }
	t_real GetBC2() const { return m_Bc2; }

	using MagSystem<t_real, t_cplx, ORDER_FOURIER>::minimise;


	std::shared_ptr<Skx<t_real, t_cplx, ORDER>> copy() const
	{
		return std::make_shared<Skx<t_real, t_cplx, ORDER>>(*this);
	}

	virtual std::shared_ptr<MagSystem<t_real, t_cplx, ORDER_FOURIER>> copyCastSys() const override { return copy(); }
	virtual std::shared_ptr<MagDynamics<t_real, t_cplx>> copyCastDyn() const override { return copy(); }


	void GenFullFourier();
	std::tuple<t_mat_cplx, t_mat_cplx> GetMCrossMFluct(
		int iGhmag, int iGkmag, t_real qh, t_real qk, t_real ql) const;
	std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
		GetSpecWeights(int iGhmag, int iGkmag, t_real qh, t_real qk, t_real ql, t_real minE=-1., t_real maxE=-2.) const;

	void SetCoords(t_real Bx, t_real By, t_real Bz,
		t_real Pinx, t_real Piny, t_real Pinz);
	void SetG(t_real h, t_real k, t_real l);
	void SetProjNeutron(bool b) { m_bProjNeutron = b; }
	void SetFilterZeroWeight(bool b) { m_filterzeroweight = b; }
	void SetWeightEps(t_real eps) { m_weighteps = eps; }
	const t_mat_cplx& GetNeutronProjOp() const { return m_projNeutron; }

	virtual std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
		GetDisp(t_real h, t_real k, t_real l, t_real minE=-1., t_real maxE=-2.) const override;


protected:
	ublas::matrix<t_vec_cplx> GetFullFourier() const;


private:
	bool m_filterzeroweight = false;
	t_real m_eveps = 1e-6;
	t_real m_weighteps = 1e-6;
	t_real m_evlimit = 0.9995;

	// B matrix with 120 deg between (100) and (010), Q_lab = B*Q_rlu
	t_mat m_Bmat = tl2::make_mat<t_mat>(
		{{1, std::cos(tl2::d2r<t_real>(120))}, {0, std::sin(tl2::d2r<t_real>(120))}});
	t_mat m_Binv;

	std::vector<t_mat> m_rot, m_rotRlu;

	t_real m_pitch = 1;
	t_real m_B = 0;
	t_real m_T = -100;
	t_real m_Bc2 = 0;

	std::vector<t_vec_cplx> m_fourier;

	// all six sixth
	std::vector<std::vector<t_vec>> m_peaks_360;
	// one sixth of the magnetic lattice
	std::vector<t_vec> m_peaks_60, m_peaks_60_lab;
	std::vector<t_vec> m_allpeaks;

	ublas::matrix<t_vec_cplx> m_M;
	t_vec m_Grlu;
	t_quat m_rotCoord;

	std::vector<t_mat_cplx> m_polMat =
	{{
		get_chiralpol<t_mat_cplx>(1),	// SF1
		get_chiralpol<t_mat_cplx>(2),	// SF2
		get_chiralpol<t_mat_cplx>(3)	// NSF
	}};

	t_mat_cplx m_projNeutron = tl2::unit_m<t_mat_cplx>(3);
	bool m_bProjNeutron = true;

	std::vector<std::pair<int, int>> m_idx1, m_idx2[3], m_idx3[4];
};


#endif
