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
 *	- The descriptions and Mathematica implementations of the different helimagnon model versions by M. Garst and J. Waizner, 2014-2018.
 *	- The 2015 and 2016 Python optimised implementations by G. Brandl and M. Kugler of the first version of the helimagnon model.
 *	  This present version started as a C++ port of G. Brandl's and M. Kugler's Python implementation,
 *	  that was then adapted to new theoretical model revisions provided by M. Garst.
 *	- The descriptions and Mathematica implementations of the different skyrmion model versions by M. Garst and J. Waizner, 2016-2020,
 *	  that also included improvements to the helimagnon model and code.
 *	- The 2016 optimised Python implementations by M. Kugler and G. Brandl of the first version of the skyrmion model,
 *	  that also included improvements to the helimagnon code.
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <cmath>
#include <numeric>
#include <functional>
#include <iostream>

#include "heli.h"
#include "heli_inst.cxx"


/**
 * constructor, precalculate loop indices
 */
template<class t_real, class t_cplx, int ORDER>
Heli<t_real, t_cplx, ORDER>::Heli()
{
	for(int i=0; i<3; ++i)
	{
		if(i<2) m_idx2[i].reserve(SIZE * SIZE);
		m_idx3[i].reserve(SIZE * SIZE * SIZE);
	}

	for(int k=-ORDER; k<=ORDER; ++k)
	{
		for(int j=-ORDER; j<=ORDER; ++j)
		{
			for(int i=-ORDER; i<=ORDER; ++i)
			{
				// unrolled indices for three loops
				int l = i - j - k;
				if(std::abs(l) <= ORDER)
				{
					m_idx3[0].push_back(i);
					m_idx3[1].push_back(j);
					m_idx3[2].push_back(k);
					m_idx3[3].push_back(l);
				}

				// unrolled indices for two loops
				int m = i - j;
				if(k == 0 && std::abs(m) <= ORDER)
				{
					m_idx2[0].push_back(i);
					m_idx2[1].push_back(j);
					m_idx2[2].push_back(m);
				}
			}
		}
	}
}


/**
 * free energy, equ. 3.37 in [waizner17]
 */
template<class t_real, class t_cplx, int ORDER>
t_real Heli<t_real, t_cplx, ORDER>::F()
{
	const t_vec_cplx& m0 = m_fourier[0];
	const auto m0_sq = tl2::inner(m0, m0);

	// dipolar interaction involving m0
	const t_mat_cplx demag = tl2::diag_matrix<t_mat_cplx>({1./3., 1./3., 1./3.});
	t_cplx cF = g_chi<t_real> * tl2::inner(m0, tl2::prod_mv(demag, m0));

	cF += (m_T_theo + 1.) * m0_sq;  // m^2 involving m0
	cF += m0_sq * m0_sq;            // m^4 involving m0

	const t_real mult = 2.;   // 2 * top peaks
	for(int i=/*-ORDER*/ 1; i<=ORDER; ++i)
	{
		const t_real q = t_real(i); // q_vec = [0, 0, q]
		const t_real q_sq = q*q;

		const t_vec_cplx& mi = get_comp(m_fourier, i);
		t_vec_cplx mj = tl2::conjugate_vec(mi);
		const auto m_sq = tl2::inner(mi, mj);

		// dipolar interaction, equ. 3.24 in [waizner17]
		cF += mult * g_chi<t_real> * mi[2]/**q*/ * mj[2]/**q / q_sq*/;

		// dmi, mi * ([0,0,q] x mj) = mi * (-q*mj[1], q*mj[0], 0])
		cF += -mult * t_cplx(0., 2.) * (-mi[0]*q*mj[1] + mi[1]*q*mj[0]);

		cF += mult * m_sq * q_sq;             // (del m)^2
		cF += mult * (m_T_theo + 1.) * m_sq;  // m^2
		cF += mult * m0_sq * m_sq;            // m^4

		// high-order correction
		cF += mult * g_hoc_b<t_real, HELI_USE_HOC> * m_sq * q_sq*q_sq;
	}

	for(std::size_t i=0; i<m_idx2[0].size(); ++i)  // m^4 involving m0
	{
		if(-m_idx2[0][i] <= 0) continue;

		const auto& m1 = get_comp(m_fourier, -m_idx2[0][i]);
		const auto& m2 = get_comp(m_fourier, +m_idx2[1][i]);
		const auto& m3 = get_comp(m_fourier, +m_idx2[2][i]);

		cF += mult * tl2::inner(m0, m1) * tl2::inner(m2, m3);
	}
	for(std::size_t i=0; i<m_idx3[0].size(); ++i)  // m^4
	{
		if(-m_idx3[0][i] <= 0) continue;

		const auto& m1 = get_comp(m_fourier, -m_idx3[0][i]);
		const auto& m2 = get_comp(m_fourier, +m_idx3[1][i]);
		const auto& m3 = get_comp(m_fourier, +m_idx3[2][i]);
		const auto& m4 = get_comp(m_fourier, +m_idx3[3][i]);

		cF += mult * tl2::inner(m1, m2) * tl2::inner(m3, m4);
	}

	cF += -m_B_theo * std::sqrt(m0_sq);  // zeeman shift
	return cF.real();
}


/**
 * set fourier components
 */
template<class t_real, class t_cplx, int ORDER>
void Heli<t_real, t_cplx, ORDER>::SetFourier(const std::vector<t_vec_cplx> &fourier, bool symm)
{
	m_fourier.reserve(SIZE);
	m_fourier = fourier;

	if(symm) // symmetrise
	{
		for(t_vec_cplx& _fourier : m_fourier)
			_fourier[1] = std::conj(_fourier[0]);
	}

	while(m_fourier.size() <= ORDER_FOURIER)
		m_fourier.emplace_back(tl2::make_vec<t_vec_cplx>({0., 0., 0.}));
	while(m_fourier.size() > ORDER_FOURIER+1)
		m_fourier.pop_back();

	// add complex conjugate bottom peaks
	for(std::size_t i=ORDER_FOURIER; i>=1; --i)
		m_fourier.push_back(tl2::conjugate_vec(m_fourier[i]));
}


/**
 * energies and spectral weights, see pp. 64-68 in [waizner17]
 */
template<class t_real, class t_cplx, int ORDER>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
Heli<t_real, t_cplx, ORDER>::GetSpecWeights(t_real qh, t_real qk, t_real ql,
	const t_mat_cplx& projNeutron, t_real minE, t_real maxE) const
{
	avoid_G(qh, qk, ql, m_eps);
	ql = -ql;

	constexpr const t_cplx imag = t_cplx(0, 1);
	t_mat_cplx Mx2d, Fluc2d;
	t_real w_scale = 1.;
	std::size_t mxrowbegin = 0;
	std::vector<t_mat_cplx> polMat(3);

	auto get_demag = [this](t_real Qi, t_real Qj, t_real Q_sq) -> t_real
	{ // equ. 6.14 in [waizner17]
		if(tl2::float_equal<t_real>(Q_sq, 0., m_eps))
			return 0.;
		return g_chi<t_real>/Q_sq * Qi*Qj;
	};

	if(!m_explicitcalc)  // calculate using ground state magnetisation
	{
		// M-cross tensor, p. 67 in [waizner17]
		auto Mx = std::make_unique<std::array<t_mat_cplx, SIZE*SIZE>>();

		for(std::size_t i=0; i<m_idx2[2].size(); ++i)
		{
			const t_vec_cplx& vecM = get_comp(m_fourier, m_idx2[2][i]);
			get_comp(*Mx, SIZE, m_idx2[0][i], m_idx2[1][i]) =
				tl2::skew<t_mat_cplx>(vecM);
		}

		// fluctuation tensor, equ. 6.13 in [waizner17]
		auto Fluc = std::make_unique<std::array<t_mat_cplx, SIZE*SIZE>>();

		for(std::size_t i=0; i<m_idx3[3].size(); ++i)
		{
			const t_vec_cplx& vecM1 = get_comp(m_fourier, m_idx3[2][i]);
			const t_vec_cplx& vecM2 = get_comp(m_fourier, m_idx3[3][i]);

			t_mat_cplx mat = 8. * tl2::outer(vecM1, vecM2) +
				4. * tl2::diag_matrix<t_mat_cplx>(3, tl2::inner(vecM1, vecM2));
			t_mat_cplx& fluccomp = get_comp(*Fluc, SIZE, m_idx3[0][i], m_idx3[1][i]);
			assign_or_add(fluccomp, mat);
		}

		for(int pk_kh=-ORDER; pk_kh<=ORDER; ++pk_kh)
		{
			t_vec Q = tl2::make_vec<t_vec>({ qh, qk, ql + t_real(pk_kh) });
			const t_real Q_sq = tl2::inner(Q, Q);

			t_mat_cplx mat(3, 3);
			for(int i=0; i<3; ++i)
			{
				// diagonal part of equ. between 6.13 and 6.14 in [waizner17]
				mat(i, i) = get_demag(Q[i], Q[i], Q_sq) + 1. + m_T_theo +
					Q_sq + g_hoc_b<t_real, HELI_USE_HOC>*Q_sq*Q_sq;

				// off-diagonal part of equ. between 6.13 and 6.14 in [waizner17]
				for(int j=i+1; j<3; ++j)
				{
					int k = 3 - i - j;               // third index in {0,1,2}
					t_real eps = (k==1 ? -1. : 1.);  // even permutation of {i,j,k}?, here: +-+
					mat(i, j) = get_demag(Q[i], Q[j], Q_sq) - eps*t_cplx(0., 2.)*Q[k];
					mat(j, i) = std::conj(mat(i, j));
				}
			}

			assign_or_add(get_comp(*Fluc, SIZE, pk_kh, pk_kh), 2.*mat);
		}

		// M-cross and fluctuation matrix
		auto mk_2dim = [](const decltype(*Mx)& arr) -> t_mat_cplx
		{
			t_mat_cplx mat = tl2::zero_matrix<t_mat_cplx>(3*SIZE, 3*SIZE);
			for(int h_idx=0; h_idx<SIZE; ++h_idx)
			{
				int h = abs_to_rel_idx(h_idx, ORDER);

				for(int k_idx=0; k_idx<SIZE; ++k_idx)
				{
					int k = abs_to_rel_idx(k_idx, ORDER);

					const t_mat_cplx& comp = get_comp(arr, SIZE, h, k);
					tl2::submatrix_copy(mat, comp, h_idx*3, k_idx*3);
				}
			}
			return mat;
		};

		Mx2d = mk_2dim(*Mx) / m_Bc2_theo;
		Fluc2d = mk_2dim(*Fluc);
		w_scale = 1.;
		mxrowbegin = 0;
		polMat[0] = get_chiralpol<t_mat_cplx>(1);   // SF1
		polMat[1] = get_chiralpol<t_mat_cplx>(2);   // SF2
		polMat[2] = get_chiralpol<t_mat_cplx>(3);   // NSF
	}
	else  // calculate using closed form
	{
		constexpr const t_real A1 = g_hoc<t_real>, A2 = A1*A1, A3 = A2*A1;
		constexpr const t_real interact = 88.;

		static const t_real Dqmin = (-2.*imag*std::pow(2., 2./3.) * std::pow(3., 5./6.) * A1 *
			std::pow(-9.*A2 + std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 1./3.) /
			(std::pow(2., 1./3.) * (3.+imag*std::sqrt(3.)) * A1 + std::pow(3., 1./6.) * (std::sqrt(3.)-imag) *
			std::pow(-9.*A2 + std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.))).real();
		static const t_real pitch_qmin = ((324.*A3 - 9.*imag*(std::sqrt(3.)-imag)*A2*std::pow(-54.*A2+6.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 1./3.) +
			(3.*imag + std::sqrt(3)) * std::sqrt(t_cplx(A3*(2.+27.*A1))) * std::pow(-54.*A2+6.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 1./3.) -
			A1*(36.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))) -3.*imag*std::pow(3., 1./6.) * std::pow(-18.*A2 + 2.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.) +
			std::pow(-54.*A2 + 6.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.))) /
			(2.*std::pow(2., 1./3.) * std::pow(3., 1./6.) * std::pow(-9.*A2+std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.) *
			(std::pow(2., 1./3.) * (-3.*imag+std::sqrt(3.))*A1 - imag*std::pow(3., 1./6) *
			(-imag+std::sqrt(3.)) * std::pow(-9.*A2+std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.)))).real();
#ifdef HELI_USE_HOC
		static const t_real hoc_qmin = (std::pow(std::pow(2., 1./3.) * (std::sqrt(3.)-3.*imag) * A1 -
			imag*std::pow(3., 1./6.)*(std::sqrt(3.)-imag) *
			std::pow(-9.*A2 + std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.), 2.) /
			(24.*std::pow(2., 1./3.) * A1*std::pow(-27.*A2+3.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.))).real();
#else
		static const t_real hoc_qmin = 0.;
#endif

		const t_real M_amp = m_B / m_Bc2;
		const t_real heli_amp = std::sqrt(0.5 - 0.5*M_amp*M_amp);
		const t_cplx hx = qh - imag*qk;

		Mx2d = tl2::zero_m<t_mat_cplx>(3*SIZE, 3*SIZE);
		Fluc2d = tl2::zero_m<t_mat_cplx>(3*SIZE, 3*SIZE);
		w_scale = g_g<t_real>;
		mxrowbegin = 3 * ORDER;

		for(int pk=0; pk<SIZE; ++pk)
		{
			const t_real Qz = ql + t_real(pk) - t_real(ORDER);
			const t_real Q_sq = qh*qh + qk*qk + Qz*Qz;

			// M-cross matrix
			Mx2d(3*pk + 0, 3*pk + 0) = -imag*M_amp;
			Mx2d(3*pk + 1, 3*pk + 1) = imag*M_amp;

			if(pk > 0)
			{
				Mx2d(3*pk + 2, 3*(pk-1) + 0) = imag*heli_amp;
				Mx2d(3*pk + 1, 3*(pk-1) + 2) = -imag*heli_amp;
			}
			if(pk < SIZE-1)
			{
				Mx2d(3*pk + 0, 3*(pk+1) + 2) = imag*heli_amp;
				Mx2d(3*pk + 2, 3*(pk+1) + 1) = -imag*heli_amp;
			}

			// fluctuation matrix
			Fluc2d(3*pk + 0, 3*pk + 0) = pitch_qmin + interact/3.*heli_amp*heli_amp
				+ 0.5*get_demag(qh*qh + qk*qk, 1., Q_sq) + 2.*Dqmin*Qz
				+ Q_sq + hoc_qmin*Q_sq*Q_sq;
			Fluc2d(3*pk + 0, 3*pk + 1) = 0.5*get_demag(1., 1., Q_sq) * hx * hx;
			Fluc2d(3*pk + 0, 3*pk + 2) = (0.5*get_demag(Qz, 1., Q_sq) - Dqmin) * std::sqrt(2) * hx;

			Fluc2d(3*pk + 1, 3*pk + 0) = std::conj(Fluc2d(3*pk + 0, 3*pk + 1));
			Fluc2d(3*pk + 1, 3*pk + 1) = Fluc2d(3*pk + 0, 3*pk + 0) - 4.*Dqmin*Qz;
			Fluc2d(3*pk + 1, 3*pk + 2) = (0.5*get_demag(Qz, 1., Q_sq) + Dqmin) * std::sqrt(2.) * std::conj(hx);

			Fluc2d(3*pk + 2, 3*pk + 0) = std::conj(Fluc2d(3*pk + 0, 3*pk + 2));
			Fluc2d(3*pk + 2, 3*pk + 1) = std::conj(Fluc2d(3*pk + 1, 3*pk + 2));
			Fluc2d(3*pk + 2, 3*pk + 2) = pitch_qmin + interact/3.*M_amp*M_amp
				+ get_demag(Qz*Qz, 1., Q_sq) + Q_sq + hoc_qmin*Q_sq*Q_sq;

			if(pk > 0)
				Fluc2d(3*pk + 1, 3*(pk-1) + 2) = Fluc2d(3*pk + 2, 3*(pk-1) + 0) = interact/3.*M_amp*heli_amp;
			if(pk < SIZE-1)
				Fluc2d(3*pk + 2, 3*(pk+1) + 1) = Fluc2d(3*pk + 0, 3*(pk+1) + 2) = interact/3.*M_amp*heli_amp;
			if(pk > 1)
				Fluc2d(3*pk + 1, 3*(pk-2) + 0) = interact/3.*heli_amp*heli_amp;
			if(pk < SIZE-2)
				Fluc2d(3*pk + 0, 3*(pk+2) + 1) = interact/3.*heli_amp*heli_amp;
		}

		polMat[0] = get_polmat<t_mat_cplx>(2); // SF1
		polMat[1] = get_polmat<t_mat_cplx>(1); // SF2
		polMat[2] = get_polmat<t_mat_cplx>(3); // NSF
	}

	// energies and weights
	return calc_weights<t_mat_cplx, t_vec_cplx, t_cplx, t_real>(
		Mx2d * g_g<t_real>, Fluc2d,
		projNeutron, polMat, g_muB<t_real> * m_Bc2, w_scale,
		minE, maxE, m_eveps, -1., m_weighteps,
		m_filterzeroweight, m_onlymode, mxrowbegin);
}


/**
 * rotates field to internal [001] convention
 */
template<class t_real, class t_cplx, int ORDER>
void Heli<t_real, t_cplx, ORDER>::SetCoords(t_real Bx, t_real By, t_real Bz, t_real /*Px*/, t_real /*Py*/, t_real /*Pz*/)
{
	t_vec B = tl2::make_vec<t_vec>({ Bx, By, Bz });
	t_quat quatB = tl2::rotation_quat(B, tl2::make_vec<t_vec>({ 0, 0, 1 }));
	m_rotCoord = quatB;
}


/**
 * dispersion
 */
template<class t_real, class t_cplx, int ORDER>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
Heli<t_real, t_cplx, ORDER>::GetDisp(t_real h, t_real k, t_real l, t_real minE, t_real maxE) const
{
	t_vec Qrlu = tl2::make_vec<t_vec>({ h, k, l });
	t_vec qrlu = Qrlu - m_Grlu;
	t_vec qkh = qrlu / g_kh_rlu<t_real>(m_T);
	qkh = tl2::quat_vec_prod(m_rotCoord, qkh);

	// orthogonal 1-|Q><Q| projector for neutron scattering
	t_mat_cplx projNeutron = tl2::unit_m<t_mat_cplx>(3);
	if(m_projNeutron)
	{
		bool bInChiralBase = m_explicitcalc;

		t_vec _Qnorm = Qrlu / tl2::veclen(Qrlu);
		t_vec_cplx Qnorm = tl2::quat_vec_prod(m_rotCoord, _Qnorm);

		if(bInChiralBase)
		{
			t_mat_cplx chiral = get_chiralbasismat<t_mat_cplx, t_vec_cplx>();
			Qnorm = tl2::prod_mv<t_vec_cplx, t_mat_cplx>(chiral, Qnorm);
		}

		projNeutron -= tl2::outer_cplx<t_vec_cplx, t_mat_cplx>(Qnorm, Qnorm);

		if(bInChiralBase)
			projNeutron = tl2::conjugate_mat(projNeutron);
	}

	return GetSpecWeights(qkh[0], qkh[1], qkh[2], projNeutron, minE, maxE);
}
