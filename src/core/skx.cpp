/**
 * free energy and dynamics in the skyrmion phase
 * @author tweber@ill.fr
 * @date jul-18
 * @desc This file implements the theoretical skyrmion model by M. Garst and J. Waizner, references:
 *	- M. Garst, J. Waizner, and D. Grundler, J. Phys. D: Appl. Phys. 50, 293002 (2017), https://doi.org/10.1088/1361-6463/aa7573
 *	- J. Waizner, PhD thesis (2017), Universität zu Köln, https://kups.ub.uni-koeln.de/7937/
 *	- Personal communications with M. Garst, 2017-2020.
 * @desc This file is based on:
 *	- The descriptions and Mathematica implementations of the different skyrmion model versions by M. Garst and J. Waizner, 2016-2020,
 *	- The 2016 optimised Python implementations by M. Kugler and G. Brandl of the first version of the skyrmion model.
 *	  This present version started as a C++ port of M. Kugler's and G. Bandl's Python implementation,
 *	  that was then adapted to new theoretical model revisions provided by M. Garst.
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <cmath>
#include <numeric>
#include <functional>
#include <iostream>

#include "skx.h"
#include "skx_inst.cxx"


/**
 * generates skyrmion satellite peaks up to a given order
 */
template<class t_vec>
static std::vector<t_vec> gen_peaks(const int ORDER)
{
	using t_val = typename t_vec::value_type;

	const int SIZE = 2*ORDER + 1;
	std::vector<t_vec> pks; pks.reserve(SIZE*SIZE);

	for(int h_idx=0; h_idx<SIZE; ++h_idx)
	{
		int h = abs_to_rel_idx(h_idx, ORDER);
		for(int k_idx=0; k_idx<SIZE; ++k_idx)
		{
			int k = abs_to_rel_idx(k_idx, ORDER);
			if(std::abs(h-k) > ORDER)
				continue;
			pks.emplace_back(tl2::make_vec<t_vec>(
				{ static_cast<t_val>(h), static_cast<t_val>(k) }));
		}
	}
	return pks;
}


/**
 * constructor, precalculate loop indices
 */
template<class t_real, class t_cplx, int ORDER>
Skx<t_real, t_cplx, ORDER>::Skx()
{
	tl2::inverse(m_Bmat, m_Binv); // xtal B matrix

	// 360 degree peaks
	m_allpeaks_rlu = gen_peaks<t_vec>(ORDER);

	// peaks in 60 degree segment
	m_peaks60rlu.reserve(ORDER_FOURIER);
	m_peaks60lab.reserve(ORDER_FOURIER);
	for(const t_vec& pk_rlu : m_allpeaks_rlu)
	{
		if(!(pk_rlu[0] > pk_rlu[1] && pk_rlu[1] >= 0))
			continue;

		t_vec pk_lab = tl2::prod_mv(m_Bmat, pk_rlu);
		pk_lab /= tl2::veclen(pk_lab);

		m_peaks60rlu.emplace_back(std::move(pk_rlu));
		m_peaks60lab.emplace_back(std::move(pk_lab));
	}

	for(int i=0; i<3; ++i)
	{
		if(i<2) m_idx2[i].reserve(m_allpeaks_rlu.size() * m_allpeaks_rlu.size());
		m_idx3[i].reserve(m_allpeaks_rlu.size() * m_allpeaks_rlu.size() * m_allpeaks_rlu.size());
	}

	for(const auto& pk_k : m_allpeaks_rlu)
	{
		for(const auto& pk_j : m_allpeaks_rlu)
		{
			// unrolled indices for two loops
			int l_h2 = lattidx(pk_j[0]) - lattidx(pk_k[0]);
			int l_k2 = lattidx(pk_j[1]) - lattidx(pk_k[1]);
			if(std::abs(l_h2) <= ORDER && std::abs(l_k2) <= ORDER
				&& std::abs(l_h2-l_k2) <= ORDER)
			{
				m_idx2[0].emplace_back(std::make_pair(lattidx(pk_j[0]), lattidx(pk_j[1])));
				m_idx2[1].emplace_back(std::make_pair(lattidx(pk_k[0]), lattidx(pk_k[1])));
				m_idx2[2].emplace_back(std::make_pair(l_h2, l_k2));
			}

			// unrolled indices for three loops
			for(const auto& pk_i : m_allpeaks_rlu)
			{
				int l_h3 = lattidx(pk_i[0]) - lattidx(pk_j[0]) - lattidx(pk_k[0]);
				int l_k3 = lattidx(pk_i[1]) - lattidx(pk_j[1]) - lattidx(pk_k[1]);
				if(std::abs(l_h3) <= ORDER && std::abs(l_k3) <= ORDER
					&& std::abs(l_h3-l_k3) <= ORDER)
				{
					m_idx3[0].emplace_back(std::make_pair(lattidx(pk_i[0]), lattidx(pk_i[1])));
					m_idx3[1].emplace_back(std::make_pair(lattidx(pk_j[0]), lattidx(pk_j[1])));
					m_idx3[2].emplace_back(std::make_pair(lattidx(pk_k[0]), lattidx(pk_k[1])));
					m_idx3[3].emplace_back(std::make_pair(l_h3, l_k3));
				}
			}
		}
	}
}


/**
 * free energy
 */
template<class t_real, class t_cplx, int ORDER>
t_real Skx<t_real, t_cplx, ORDER>::F()
{
	auto is_peak_in_top_half = [this](const t_vec& q, t_real q_sq) -> bool
	{
		bool in_top = ((tl2::float_equal<t_real>(q[1], 0., m_eps) && q[0] >= 0.) || q[1] > 0.);
		bool valid_pks = !tl2::float_equal<t_real>(q_sq, 0., m_eps);

		return in_top && valid_pks;
	};

	auto is_hk_in_top_half = [this, &is_peak_in_top_half](int h, int k) -> bool
	{
		t_vec q_rlu = tl2::make_vec<t_vec>({ t_real(h), t_real(k) });
		t_vec q = tl2::prod_mv(m_Bmat, q_rlu);
		const t_real q_sq = tl2::inner(q, q);

		return is_peak_in_top_half(q, q_sq);
	};

	const t_vec_cplx& m0 = m_M(0, 0);
	const auto m0_sq = tl2::inner(m0, m0);

	// dipolar interaction
	static const t_mat_cplx demag = tl2::diag_matrix<t_mat_cplx>({1./3., 1./3., 1./3.});
	t_cplx cF = g_chi<t_real> * tl2::inner(m0, tl2::prod_mv(demag, m0));

	cF += (m_T + 1.) * m0_sq; // phi^2
	cF += m0_sq * m0_sq;      // phi^4

	const t_real mult = 2.;   // 2 * top 180 degree peaks
	for(const t_vec& q_rlu : m_allpeaks_rlu)
	{
		const int hk[2] = { lattidx(q_rlu[0]), lattidx(q_rlu[1]) };

		t_vec q = tl2::prod_mv(m_Bmat, q_rlu);
		const t_real q_sq = tl2::inner(q, q);

		if(!is_peak_in_top_half(q, q_sq))
			continue;

		t_vec_cplx qc = tl2::make_vec<t_vec_cplx>({ q[0], q[1], 0. });
		const t_vec_cplx& m = get_comp(m_M, hk[0], hk[1]);
		t_vec_cplx mj = tl2::conjugate_vec(m);
		const t_cplx m_sq = tl2::inner(m, mj);

		// dipolar interaction
		if(!tl2::float_equal<t_real>(q_sq, 0., m_eps))
			cF += mult * g_chi<t_real> * tl2::inner(m, qc) * tl2::inner(mj, qc) / q_sq;

		cF += -2. * mult * t_cplx(0, 1) * tl2::inner(m, tl2::cross_3(qc, mj)); // dmi
		cF += mult * m_sq * q_sq;        // phi^2
		cF += mult * (m_T + 1.) * m_sq;  // phi^2
		cF += mult * m0_sq * m_sq;       // phi^4
		cF += mult * g_hoc_b<t_real, SKX_USE_HOC> * m_sq * q_sq*q_sq; // high-order correction
	}

	for(std::size_t i=0; i<m_idx2[0].size(); ++i)  // phi^4
	{
		const int hk[2] = { -m_idx2[0][i].first, -m_idx2[0][i].second };
		if(!is_hk_in_top_half(hk[0], hk[1]))
			continue;

		const auto& m1 = get_comp(m_M, hk[0], hk[1]);
		const auto& m2 = get_comp(m_M, m_idx2[1][i].first, m_idx2[1][i].second);
		const auto& m3 = get_comp(m_M, m_idx2[2][i].first, m_idx2[2][i].second);

		cF += mult * tl2::inner(m0, m1) * tl2::inner(m2, m3);
	}
	for(std::size_t i=0; i<m_idx3[0].size(); ++i)  // phi^4
	{
		const int hk[2] = { -m_idx3[0][i].first, -m_idx3[0][i].second };
		if(!is_hk_in_top_half(hk[0], hk[1]))
			continue;

		const auto& m1 = get_comp(m_M, hk[0], hk[1]);
		const auto& m2 = get_comp(m_M, m_idx3[1][i].first, m_idx3[1][i].second);
		const auto& m3 = get_comp(m_M, m_idx3[2][i].first, m_idx3[2][i].second);
		const auto& m4 = get_comp(m_M, m_idx3[3][i].first, m_idx3[3][i].second);

		cF += mult * tl2::inner(m1, m2) * tl2::inner(m3, m4);
	}

	cF += -m_B * std::sqrt(m0_sq);  // zeeman shift
	return cF.real();
}


/**
 * set fourier components
 */
template<class t_real, class t_cplx, int ORDER>
void Skx<t_real, t_cplx, ORDER>::SetFourier(const std::vector<t_vec_cplx> &fourier, bool symm)
{
	m_fourier.reserve(ORDER_FOURIER+1);
	m_fourier = fourier;
	if(symm) // symmetrise
	{
		for(t_vec_cplx& _fourier : m_fourier)
			_fourier[1] = std::conj(_fourier[0]);
	}

	while(m_fourier.size() < ORDER_FOURIER+1)
		m_fourier.emplace_back(tl2::make_vec<t_vec_cplx>({0., 0., 0.}));
	while(m_fourier.size() > ORDER_FOURIER+1)
		m_fourier.pop_back();

	// generate full fourier coefficients
	m_M = tl2::make_mat<ublas::matrix<t_vec_cplx>>(SIZE, SIZE, tl2::make_vec<t_vec_cplx>({0, 0, 0}));
	m_M(0,0) = m_fourier[0];

	// generate all skx fourier components
	for(int rot_idx=0; rot_idx<6; ++rot_idx)
	{
		const t_mat rotLab = tl2::rotation_matrix_2d<t_mat>(tl2::d2r<t_real>(60*rot_idx));
		const t_mat rotRlu = tl2::transform<t_mat>(rotLab, m_Bmat, 0);

		for(std::size_t peak_idx=0; peak_idx<m_peaks60rlu.size(); ++peak_idx)
		{
			const t_vec pk_rlu = tl2::prod_mv(rotRlu, m_peaks60rlu[peak_idx]);
			const int hk[2] = { lattidx(pk_rlu[0]), lattidx(pk_rlu[1]) };

			t_vec_cplx& M = get_comp(m_M, hk[0], hk[1]);
			t_vec_cplx m(2);

			if(hk[0] != 0 || hk[1] != 0)  // avoid singularity at (0, 0)
			{
				m[0] = m_peaks60lab[peak_idx][1] * m_fourier[peak_idx+1][1];
				m[1] = m_peaks60lab[peak_idx][0] * m_fourier[peak_idx+1][0];
				m = tl2::prod_mv(rotLab, m);
			}

			M[0] = m[0]; M[1] = m[1]; M[2] = m_fourier[peak_idx+1][2];
		}
	}
}


/**
 * energies and spectral weights
 */
template<class t_real, class t_cplx, int ORDER>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
Skx<t_real, t_cplx, ORDER>::GetSpecWeights(int Ghmag, int Gkmag,
	t_real qh, t_real qk, t_real ql, const t_mat_cplx& projNeutron,
	t_real minE, t_real maxE) const
{
	avoid_G(qh, qk, ql, m_eps);
	const t_vec q_lab = tl2::make_vec<t_vec>({ qh, qk, -ql });

	// M-cross tensor
	auto Mx = std::make_unique<std::array<t_mat_cplx, SIZE*SIZE*SIZE*SIZE>>();

	for(std::size_t i=0; i<m_idx2[2].size(); ++i)
	{
		get_comp(*Mx, SIZE,
			m_idx2[0][i].first, m_idx2[0][i].second,
			m_idx2[1][i].first, m_idx2[1][i].second) =
				tl2::skew<t_mat_cplx>(
					get_comp(m_M, m_idx2[2][i].first, m_idx2[2][i].second));
	}

	// fluctuation tensor
	auto Fluc = std::make_unique<std::array<t_mat_cplx, SIZE*SIZE*SIZE*SIZE>>();

	for(std::size_t i=0; i<m_idx3[3].size(); ++i)
	{
		const t_vec_cplx& vecM1 = get_comp(m_M, m_idx3[2][i].first, m_idx3[2][i].second);
		const t_vec_cplx& vecM2 = get_comp(m_M, m_idx3[3][i].first, m_idx3[3][i].second);

		t_mat_cplx mat = 8.*tl2::outer(vecM1, vecM2) +
			4.*tl2::diag_matrix<t_mat_cplx>(3, tl2::inner(vecM1, vecM2));
		t_mat_cplx& fluccomp = get_comp(*Fluc, SIZE,
			m_idx3[0][i].first, m_idx3[0][i].second,
			m_idx3[1][i].first, m_idx3[1][i].second);
		assign_or_add(fluccomp, mat);
	}

	for(const t_vec &pk_rlu : m_allpeaks_rlu)
	{
		t_vec pk_lab = tl2::prod_mv(m_Bmat, pk_rlu);
		pk_lab.resize(3, true); pk_lab[2] = 0.;
		const t_vec Q = pk_lab + q_lab;
		const t_real Q_sq = tl2::inner(Q, Q);

		auto get_dip = [this](t_real Qi, t_real Qj, t_real Q_sq) -> t_real
		{
			if(tl2::float_equal<t_real>(Q_sq, 0., m_eps))
				return 0.;
			return g_chi<t_real>/Q_sq * Qi*Qj;
		};

		t_mat_cplx mat(3,3);
		for(int i=0; i<3; ++i)  // diagonal
			mat(i, i) = get_dip(Q[i], Q[i], Q_sq) + 1. + m_T + Q_sq
				+ g_hoc_b<t_real, SKX_USE_HOC>*Q_sq*Q_sq;
		for(int i=0; i<2; ++i)  // off-diagonal
		{
			for(int j=i+1; j<3; ++j)
			{
				int k = 3 - i - j;                // third index in {0,1,2}
				t_real sign = (k==1 ? 1. : -1.);  // - + -
				mat(i, j) = get_dip(Q[i], Q[j], Q_sq) + sign*2.*t_cplx(0,1)*Q[k];
				mat(j, i) = std::conj(mat(i, j));
			}
		}

		const int hk[2] = { lattidx(pk_rlu[0]), lattidx(pk_rlu[1]) };
		assign_or_add(get_comp(*Fluc, SIZE, hk[0], hk[1], hk[0], hk[1]), 2.*mat);
	}

	// M-cross and fluctuation matrix
	const int MAXORDER = ORDER + std::max(std::abs(Ghmag), std::abs(Gkmag));

	auto mk_2dim = [MAXORDER, Ghmag, Gkmag](const decltype(*Mx)& arr) -> t_mat_cplx
	{
		std::vector<t_vec_int> pks = gen_peaks<t_vec_int>(MAXORDER);

		t_mat_cplx mat = tl2::zero_matrix<t_mat_cplx>(3*pks.size(), 3*pks.size());
		for(std::size_t h_idx=0; h_idx<pks.size(); ++h_idx)
		{
			for(std::size_t k_idx=0; k_idx<pks.size(); ++k_idx)
			{
				const t_mat_cplx& comp = get_flat_comp(
					arr, SIZE, ORDER, MAXORDER,
					pks[h_idx][0] + Ghmag, pks[h_idx][1] + Gkmag,
					pks[k_idx][0] + Ghmag, pks[k_idx][1] + Gkmag);

				tl2::submatrix_copy(mat, comp, h_idx*3, k_idx*3);
			}
		}
		return mat;
	};

	// energies and weights
	return calc_weights<t_mat_cplx, t_vec_cplx, t_cplx, t_real>(
		mk_2dim(*Mx) * g_g<t_real> / m_Bc2, mk_2dim(*Fluc),
		projNeutron, m_polMat, g_muB<t_real> * m_Bc2_exp, 1.,
		minE, maxE, m_eveps, m_evlimit, m_weighteps,
		m_filterzeroweight, /*m_onlymode*/-1);
}


/**
 * rotates field and pinning to internal conventions ([001] and [100], respectively)
 */
template<class t_real, class t_cplx, int ORDER>
void Skx<t_real, t_cplx, ORDER>::SetCoords(
	t_real Bx, t_real By, t_real Bz,
	t_real Pinx, t_real Piny, t_real Pinz)
{
	t_vec B = tl2::make_vec<t_vec>({ Bx, By, Bz });
	t_vec _Pin = tl2::make_vec<t_vec>({ Pinx, Piny, Pinz });

	t_quat quatB = tl2::rotation_quat(B, tl2::make_vec<t_vec>({ 0, 0, 1 }));

	t_vec Pin = tl2::quat_vec_prod(quatB, _Pin);
	t_quat quatPin = tl2::rotation_quat(Pin, tl2::make_vec<t_vec>({ 1, 0, 0 }));

	m_rotCoord = quatPin * quatB;
}


/**
 * query the dispersion
 */
template<class t_real, class t_cplx, int ORDER>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
Skx<t_real, t_cplx, ORDER>::GetDisp(t_real h, t_real k, t_real l, t_real minE, t_real maxE) const
{
	t_vec Qrlu = tl2::make_vec<t_vec>({ h, k, l });
	t_vec qrlu = Qrlu - m_Grlu;
	t_vec qkh = qrlu / g_kh_rlu_29K<t_real>;

	// orthogonal 1-|Q><Q| projector for neutron scattering
	t_mat_cplx projNeutron = tl2::unit_m<t_mat_cplx>(3);
	if(m_bProjNeutron)
	{
		t_vec Qnorm = Qrlu / tl2::veclen(Qrlu);
		Qnorm = tl2::quat_vec_prod(m_rotCoord, Qnorm);
		projNeutron -= tl2::outer<t_vec, t_mat>(Qnorm, Qnorm);
	}

	qkh = tl2::quat_vec_prod(m_rotCoord, qkh);
	t_real _l = qkh[2];
	qkh.resize(2, true);

	t_vec Qmagrlu = tl2::prod_mv(m_Binv, qkh);
	t_vec Gmagrlu = tl2::make_vec<t_vec>({ std::round(Qmagrlu[0]), std::round(Qmagrlu[1]) });

	const std::vector<t_vec> sats = gen_peaks<t_vec>(1);

	auto iterClosest = std::min_element(sats.begin(), sats.end(),
		[&Gmagrlu, &Qmagrlu, this](const t_vec& sat1, const t_vec& sat2) -> bool
		{
			t_vec qmag1 = Qmagrlu - (Gmagrlu+sat1);
			t_vec qmag2 = Qmagrlu - (Gmagrlu+sat2);

			t_vec q1 = tl2::prod_mv(m_Bmat, qmag1);
			t_vec q2 = tl2::prod_mv(m_Bmat, qmag2);

			return tl2::veclen_sq(q1) < tl2::veclen_sq(q2);
		});

	Gmagrlu += *iterClosest;
	t_vec qmagrlu = Qmagrlu - Gmagrlu;
	t_vec qmaglab = tl2::prod_mv(m_Bmat, qmagrlu);

	return GetSpecWeights(int(Gmagrlu[0]), int(Gmagrlu[1]),
		qmaglab[0], qmaglab[1], _l, projNeutron, minE, maxE);
}
