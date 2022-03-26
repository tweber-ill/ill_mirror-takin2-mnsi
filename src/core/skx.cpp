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
 *	  This present version started as a C++ port of that Python implementation by M. Kugler and G. Brandl,
 *	  that was then adapted to new theoretical model revisions provided by M. Garst.
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <cmath>
#include <numeric>
#include <functional>
#include <iostream>

#include "skx.h"


// instantiation
#ifdef DEF_SKX_ORDER
	#pragma message("Skx Order: " __TL2_STRCONV(DEF_SKX_ORDER))
	template class Skx<double, std::complex<double>, DEF_SKX_ORDER>;

	#ifdef __HACK_FULL_INST__
		template Skx<double, std::complex<double>, DEF_SKX_ORDER>::Skx();
		template void Skx<double, std::complex<double>, DEF_SKX_ORDER>::SetG(double, double, double);
		template void Skx<double, std::complex<double>, DEF_SKX_ORDER>::SetCoords(
			double, double, double, double, double, double);
	#endif
#endif


/**
 * constructor
 */
template<class t_real, class t_cplx, int ORDER>
Skx<t_real, t_cplx, ORDER>::Skx()
{
	// xtal B matrix
	tl2::inverse(m_Bmat, m_Binv);

	// rotation in lab
	t_mat curRot = tl2::unit_m<t_mat>(2);
	const auto rot60 = tl2::rotation_matrix_2d<t_mat>(tl2::d2r<t_real>(60));

	m_rot.reserve(6);
	for(int i=0; i<6; ++i)
	{
		curRot = tl2::prod_mm(curRot, rot60);
		m_rot.push_back(curRot);
	}

	// rotation in rlu
	m_rotRlu.reserve(m_rot.size());
	for(const auto& rot : m_rot)
	{
		auto M = tl2::transform<t_mat>(rot, m_Bmat, 0);
		m_rotRlu.emplace_back(M);
	}


	// all peaks
	m_allpeaks.reserve((2*ORDER+1) * (2*ORDER+1));
	for(t_real h=-ORDER; h<ORDER+1; ++h)
	{
		for(t_real k=-ORDER; k<ORDER+1; ++k)
		{
			if(std::abs(h-k) <= t_real(ORDER))
				m_allpeaks.emplace_back(tl2::make_vec<t_vec>({h,k}));
		}
	}

	m_peaks_60.reserve((ORDER+1) * ORDER);
	m_peaks_60_lab.reserve((ORDER+1) * ORDER);
	for(t_real h=0; h<ORDER+1; ++h)
	{
		for(t_real k=0; k<h; ++k)
		{
			t_vec vec_rlu = tl2::make_vec<t_vec>({h, k});
			t_vec vec_lab = tl2::prod_mv(m_Bmat, vec_rlu);
			vec_lab /= tl2::veclen(vec_lab);

			m_peaks_60.emplace_back(std::move(vec_rlu));
			m_peaks_60_lab.emplace_back(std::move(vec_lab));
		}
	}


	// indices for the top 180 degrees
	m_idx_top.reserve(m_rotRlu.size() * m_peaks_60.size());

	m_peaks_360.reserve(m_rotRlu.size());
	for(const auto& rot : m_rotRlu)
	{
		std::vector<t_vec> _peaks;
		_peaks.reserve(m_peaks_60.size());

		for(const auto &vec : m_peaks_60)
		{
			t_vec pk_rlu = tl2::prod_mv(rot, vec);
			t_vec pk = tl2::prod_mv(m_Bmat, pk_rlu);

			if((tl2::float_equal<t_real>(pk[1], 0.) && pk[0] >= 0.)
				|| pk[1] > 0.)
			{
				m_idx_top.emplace_back(std::make_pair(
					int(std::round(pk_rlu[0])),
					int(std::round(pk_rlu[1]))));
			}

			_peaks.emplace_back(pk_rlu);
		}

		m_peaks_360.emplace_back(_peaks);
	}


	auto reduce_2 = [](decltype(m_idx2)& _idx2, int sign)
	{
		decltype(m_idx2) idx2;
		for(int i=0; i<3; ++i)
			idx2[i].reserve(_idx2[0].size());
		for(std::size_t i=0; i<_idx2[0].size(); ++i)
		{
			auto val1 = sign*_idx2[0][i].first - _idx2[1][i].first;
			auto val2 = sign*_idx2[0][i].second - _idx2[1][i].second;
			if(std::abs(val1) > ORDER || std::abs(val2) > ORDER || std::abs(val1-val2) > ORDER)
				continue;

			idx2[0].push_back(_idx2[0][i]);
			idx2[1].push_back(_idx2[1][i]);
			idx2[2].emplace_back(std::make_pair(val1, val2));
		}
		for(int i=0; i<3; ++i)
			_idx2[i] = std::move(idx2[i]);
	};

	auto reduce_3 = [](decltype(m_idx3)& _idx3, int sign)
	{
		decltype(m_idx3) idx3;
		for(int i=0; i<4; ++i)
			idx3[i].reserve(_idx3[0].size());
		for(std::size_t i=0; i<_idx3[0].size(); ++i)
		{
			auto val1 = sign*_idx3[0][i].first - _idx3[1][i].first - _idx3[2][i].first;
			auto val2 = sign*_idx3[0][i].second - _idx3[1][i].second - _idx3[2][i].second;
			if(std::abs(val1) > ORDER || std::abs(val2) > ORDER || std::abs(val1-val2) > ORDER)
				continue;

			idx3[0].push_back(_idx3[0][i]);
			idx3[1].push_back(_idx3[1][i]);
			idx3[2].push_back(_idx3[2][i]);
			idx3[3].emplace_back(std::make_pair(val1, val2));
		}
		for(int i=0; i<4; ++i)
			_idx3[i] = std::move(idx3[i]);
	};


	// #2
	for(int i=0; i<2; ++i)
	{
		m_idx2[i].reserve(m_allpeaks.size() * m_idx_top.size());
		m_idx2_dyn[i].reserve(m_allpeaks.size() * m_allpeaks.size());
	}

	for(std::size_t i=0; i<m_allpeaks.size(); ++i)
	{
		std::copy(m_idx_top.begin(), m_idx_top.end(), std::back_inserter(m_idx2[0]));
		for(std::size_t j=0; j<m_idx_top.size(); ++j)
		{
			m_idx2[1].emplace_back(std::make_pair(
				int(std::round(m_allpeaks[i][0])),
				int(std::round(m_allpeaks[i][1]))));
		}

		for(std::size_t j=0; j<m_allpeaks.size(); ++j)
		{
			m_idx2_dyn[0].emplace_back(std::make_pair(
				int(std::round(m_allpeaks[j][0])),
				int(std::round(m_allpeaks[j][1]))));
			m_idx2_dyn[1].emplace_back(std::make_pair(
				int(std::round(m_allpeaks[i][0])),
				int(std::round(m_allpeaks[i][1]))));
		}
	}

	// #3
	for(int i=0; i<3; ++i)
	{
		m_idx3[i].reserve(m_allpeaks.size() * m_allpeaks.size() * m_idx_top.size());
		m_idx3_dyn[i].reserve(m_allpeaks.size() * m_allpeaks.size() * m_allpeaks.size());
	}

	for(std::size_t i=0; i<m_allpeaks.size(); ++i)
	{
		std::copy(m_idx2[0].begin(), m_idx2[0].end(), std::back_inserter(m_idx3[0]));
		std::copy(m_idx2[1].begin(), m_idx2[1].end(), std::back_inserter(m_idx3[1]));

		std::copy(m_idx2_dyn[0].begin(), m_idx2_dyn[0].end(), std::back_inserter(m_idx3_dyn[0]));
		std::copy(m_idx2_dyn[1].begin(), m_idx2_dyn[1].end(), std::back_inserter(m_idx3_dyn[1]));

		for(std::size_t j=0; j<m_idx_top.size(); ++j)
		{
			for(std::size_t k=0; k<m_allpeaks.size(); ++k)
			{
				m_idx3[2].emplace_back(std::make_pair(
					int(std::round(m_allpeaks[i][0])),
					int(std::round(m_allpeaks[i][1]))));
			}
		}

		for(std::size_t j=0; j<m_allpeaks.size(); ++j)
		{
			for(std::size_t k=0; k<m_allpeaks.size(); ++k)
			{
				m_idx3_dyn[2].emplace_back(std::make_pair(
					int(std::round(m_allpeaks[i][0])),
					int(std::round(m_allpeaks[i][1]))));
			}
		}
	}

	reduce_2(m_idx2, -1);
	reduce_3(m_idx3, -1);

	reduce_2(m_idx2_dyn, 1);
	reduce_3(m_idx3_dyn, 1);
}


template<class t_real, class t_cplx, int ORDER>
ublas::matrix<typename Skx<t_real, t_cplx, ORDER>::t_vec_cplx>
Skx<t_real, t_cplx, ORDER>::GetFullFourier() const
{
	constexpr auto imag = t_cplx(0, 1);

	ublas::matrix<t_vec_cplx> M(2*ORDER+1, 2*ORDER+1);
	for(std::size_t i=0; i<M.size1(); ++i)
		for(std::size_t j=0; j<M.size2(); ++j)
			M(i,j) = tl2::make_vec<t_vec_cplx>({0,0,0});
	M(0,0) = m_fourier[0];

	// generate all skx fourier components
	for(std::size_t ipk=0; ipk<m_peaks_360.size(); ++ipk)	// 6
	{
		const auto& rot = m_rot[ipk];
		for(std::size_t ihx=0; ihx<m_peaks_360[ipk].size(); ++ihx)	// 10
		{
			int idx1 = int(std::round(m_peaks_360[ipk][ihx][0]));
			int idx2 = int(std::round(m_peaks_360[ipk][ihx][1]));

			const auto& vecPk = m_peaks_60_lab[ihx];
			auto fourier = tl2::make_vec<t_vec_cplx>
				({ -vecPk[1] * m_fourier[ihx+1][0],
				   +vecPk[0] * m_fourier[ihx+1][0] });
			fourier = tl2::prod_mv(rot, fourier);

			get_comp(M, idx1, idx2) = tl2::make_vec<t_vec_cplx>
				({ imag * fourier[0],
				   imag * fourier[1],
				   m_fourier[ihx+1][2] });
		}
	}

	return M;
}


template<class t_real, class t_cplx, int ORDER>
void Skx<t_real, t_cplx, ORDER>::GenFullFourier()
{
	m_M = GetFullFourier();
}


/**
 * free energy
 */
template<class t_real, class t_cplx, int ORDER>
t_real Skx<t_real, t_cplx, ORDER>::F()
{
	constexpr auto imag = t_cplx(0,1);

	GenFullFourier();

	// free energy
	t_cplx cF = 0;
	auto m0 = tl2::veclen(m_M(0,0));

	// dip
	cF += g_chi<t_real>/3. * m0*m0;

	for(const auto& pair : m_idx_top)
	{
		t_vec q_rlu = tl2::make_vec<t_vec>({t_real(pair.first), t_real(pair.second)});
		t_vec q = tl2::prod_mv(m_Bmat, q_rlu);

		const auto& m = get_comp(m_M, pair.first, pair.second);
		const auto m_sq = tl2::inner_cplx(m, m);
		const auto q_sq = tl2::inner(q, q);

		// dip
		cF += 2. * g_chi<t_real> * std::pow(q[0]*m[0] + q[1]*m[1], 2.) / q_sq;

		// dmi
		cF += 8. * imag * m[2] * (q[0]*m[1] - q[1]*m[0]);

		// hoc
		cF += 2. * g_hoc<t_real> * m_sq * q_sq*q_sq;

		// phi^2 & phi^4
		cF += 2. * m_sq * q_sq;
		cF += 2. * (m_T + 1. + m0*m0) * m_sq;
	}

	// phi^2 & phi^4
	cF += (m_T + 1.) * m0*m0 + m0*m0*m0*m0;

	for(std::size_t i=0; i<m_idx2[0].size(); ++i)
	{
		const auto& m1 = get_comp(m_M, m_idx2[0][i].first, m_idx2[0][i].second);
		const auto& m2 = get_comp(m_M, m_idx2[1][i].first, m_idx2[1][i].second);
		const auto& m3 = get_comp(m_M, m_idx2[2][i].first, m_idx2[2][i].second);

		// the x and y components are purely imaginary, z purely real
		cF += 2. * m0 * m1[2] * tl2::inner(m2, m3);
	}

	// phi^4
	for(std::size_t i=0; i<m_idx3[0].size(); ++i)
	{
		const auto& m1 = get_comp(m_M, m_idx3[0][i].first, m_idx3[0][i].second);
		const auto& m2 = get_comp(m_M, m_idx3[1][i].first, m_idx3[1][i].second);
		const auto& m3 = get_comp(m_M, m_idx3[2][i].first, m_idx3[2][i].second);
		const auto& m4 = get_comp(m_M, m_idx3[3][i].first, m_idx3[3][i].second);

		// the x and y components are purely imaginary, z purely real
		cF += 2. * tl2::inner(m1, m2) * tl2::inner(m3, m4);
	}

	// zee
	cF += -m_B*m0;
	return cF.real();
}


/**
 * set fourier components
 */
template<class t_real, class t_cplx, int ORDER>
void Skx<t_real, t_cplx, ORDER>::SetFourier(const std::vector<t_vec_cplx> &fourier)
{
	m_fourier = fourier;
	while(m_fourier.size() < ORDER_FOURIER+1)
		m_fourier.emplace_back(tl2::make_vec<t_vec_cplx>({0., 0., 0.}));
	while(m_fourier.size() > ORDER_FOURIER+1)
		m_fourier.pop_back();
}


/**
 * cross-product and fluctuation matrices
 */
template<class t_real, class t_cplx, int ORDER>
std::tuple<typename Skx<t_real, t_cplx, ORDER>::t_mat_cplx, typename Skx<t_real, t_cplx, ORDER>::t_mat_cplx>
Skx<t_real, t_cplx, ORDER>::GetMCrossMFluct(
	int iGh, int iGk, t_real qh, t_real qk, t_real ql) const
{
	constexpr int SIZE = 2*ORDER+1;
	constexpr auto imag = t_cplx(0,1);

	const int CRYSSIZE = ORDER + std::max(std::abs(iGh), std::abs(iGk));
	const int VIRTSIZE = 2*CRYSSIZE + 1;

	t_vec q_lab = tl2::prod_mv(m_Bmat, tl2::make_vec<t_vec>({ qh, qk }));


	// ------------------------------------------------------------------------
	// Mx_ij = eps_ikj M_k
	auto Mx = std::make_unique<std::array<t_mat_cplx, SIZE*SIZE*SIZE*SIZE>>();

	for(std::size_t i=0; i<m_idx2_dyn[2].size(); ++i)
	{
		const auto& vecM = get_comp(m_M, m_idx2_dyn[2][i].first, m_idx2_dyn[2][i].second);
		auto skew = tl2::skew<t_mat_cplx>(vecM);

		get_comp(*Mx, SIZE,
			m_idx2_dyn[0][i].first, m_idx2_dyn[0][i].second,
			m_idx2_dyn[1][i].first, m_idx2_dyn[1][i].second) = skew;
	}
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// fluct. F_q1q2_ij = 0.5 * d^2F/(dM_-q_1_q dM_q2_j)
	auto Fluc = std::make_unique<std::array<t_mat_cplx, SIZE*SIZE*SIZE*SIZE>>();

	for(std::size_t i=0; i<m_idx3_dyn[3].size(); ++i)
	{
		const auto& vecM1 = get_comp(m_M, m_idx3_dyn[2][i].first, m_idx3_dyn[2][i].second);
		const auto& vecM2 = get_comp(m_M, m_idx3_dyn[3][i].first, m_idx3_dyn[3][i].second);

		t_mat_cplx mat = 4. * tl2::outer(vecM1, vecM2);

		// delta term
		t_cplx m1m2 = tl2::inner(vecM1, vecM2);
		for(int d=0; d<3; ++d)
			mat(d,d) += 2.*m1m2;

		auto& oldmat = get_comp(*Fluc, SIZE,
			m_idx3_dyn[0][i].first, m_idx3_dyn[0][i].second,
			m_idx3_dyn[1][i].first, m_idx3_dyn[1][i].second);
		if(!oldmat.size1())
			oldmat = 2.*mat;
		else
			oldmat += 2.*mat;
	}

	for(const auto &pk_rlu : m_allpeaks)
	{
		t_vec pk_lab = tl2::prod_mv(m_Bmat, pk_rlu);

		t_vec Q = pk_lab + q_lab;
		Q.resize(3, true);
		Q[2] = ql;

		t_real Q_sq = tl2::inner(Q, Q);
		t_real dipole = g_chi<t_real> / Q_sq;

		// diagonal
		t_mat_cplx mat(3,3);
		for(int i=0; i<3; ++i)
			mat(i,i) = dipole * Q[i]*Q[i] + 1. + m_T + Q_sq;

		// upper right triangle
		for(int i=0; i<2; ++i)
		{
			for(int j=i+1; j<3; ++j)
			{
				int k = 3 - i - j; // third index in {0,1,2}
				t_real sign = (k==1 ? 1. : -1.);
				mat(i,j) = dipole * Q[i]*Q[j] + sign*2.*imag*Q[k];
				mat(j,i) = std::conj(mat(i,j));
			}
		}

		int idx1 = int(std::round(pk_rlu[0]));
		int idx2 = int(std::round(pk_rlu[1]));

		auto& oldmat = get_comp(*Fluc, SIZE, idx1, idx2, idx1, idx2);
		if(!oldmat.size1())
			oldmat = 2.*mat;
		else
			oldmat += 2.*mat;
	}
	// ------------------------------------------------------------------------


	auto mk_2dim = [VIRTSIZE, CRYSSIZE, iGh, iGk](
		const decltype(*Mx)& arr) -> t_mat_cplx
	{
		std::vector<int> pks1(VIRTSIZE);
		std::iota(pks1.begin(), pks1.begin()+CRYSSIZE+1, 0);
		std::iota(pks1.begin()+CRYSSIZE+1, pks1.end(), -CRYSSIZE);

		std::vector<std::pair<int, int>> pks2;
		pks2.reserve(VIRTSIZE*VIRTSIZE);
		for(int k : pks1)
		{
			for(int h : pks1)
			{
				if(std::abs(h-k) > CRYSSIZE)
					continue;
				pks2.emplace_back(std::make_pair(h, k));
			}
		}

		auto mat = tl2::zero_matrix<t_mat_cplx>(3*pks2.size(), 3*pks2.size());
		for(std::size_t idx1=0; idx1<pks2.size(); ++idx1)
		{
			const auto& hk1 = pks2[idx1];
			for(std::size_t idx2=0; idx2<pks2.size(); ++idx2)
			{
				const auto& hk2 = pks2[idx2];
				const auto& comp = get_virt_comp(
					arr, SIZE, VIRTSIZE, ORDER,
					hk1.first+iGh, hk1.second+iGk,
					hk2.first+iGh, hk2.second+iGk);

				tl2::submatrix_copy(mat, comp, idx1*3, idx2*3);
			}
		}

		return mat;
	};


	auto Mx2d = mk_2dim(*Mx);
	auto Fluc2d = mk_2dim(*Fluc);
	return std::make_tuple(Mx2d, Fluc2d);
}


/**
 * energies and spectral weights
 */
template<class t_real, class t_cplx, int ORDER>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
Skx<t_real, t_cplx, ORDER>::GetSpecWeights(
	int iGhmag, int iGkmag, t_real qh, t_real qk, t_real ql,
	t_real minE, t_real maxE) const
{
	const t_real epshkl = 1e-5;
	if(tl2::float_equal<t_real>(qh, 0., epshkl) &&
		tl2::float_equal<t_real>(qk, 0., epshkl) &&
		tl2::float_equal<t_real>(ql, 0., epshkl))
	{
		qh += epshkl;
		qk += epshkl;
		ql += epshkl;
	}

	ql = -ql;

	t_real sqrtfac = std::sqrt(t_real(-0.5 - m_T*0.5));
	// energy scaling depends on Hc2_int (and the sample's demagnetisation factor)
	// Hc2_int changes rapidly vs. T here at the border to the paramagnetic phase
	const t_real E_scale_fac_heli = 0.0387;		// calculated with heli.cpp
	t_real E_scale_fac = E_scale_fac_heli / (sqrtfac/10.);

	// set experimental value for Hc2_int if given
	if(m_Bc2_exp >= 0.)
		E_scale_fac = g_muB<t_real> * m_Bc2_exp;

	t_mat_cplx Mx2d, Fluc2d;
	std::tie(Mx2d, Fluc2d) = GetMCrossMFluct(iGhmag, iGkmag, qh, qk, ql);

	// energies and weights
	return calc_weights<t_mat_cplx, t_vec_cplx, t_cplx, t_real>(
		Mx2d, Fluc2d,
		m_bProjNeutron, m_projNeutron, m_polMat,
		sqrtfac, E_scale_fac,
		minE, maxE,
		m_eveps, m_evlimit, m_weighteps,
		m_filterzeroweight, /*m_onlymode*/-1);
}


/**
 * rotates field and pinning to internal conventions
 * ([001] and [100], respectively])
 */
template<class t_real, class t_cplx, int ORDER>
void Skx<t_real, t_cplx, ORDER>::SetCoords(
	t_real Bx, t_real By, t_real Bz,
	t_real Pinx, t_real Piny, t_real Pinz)
{
	t_vec B = tl2::make_vec<t_vec>( {Bx, By, Bz} );
	t_vec _Pin = tl2::make_vec<t_vec>( {Pinx, Piny, Pinz} );

	t_quat quatB = tl2::rotation_quat(B, tl2::make_vec<t_vec>( {0, 0, 1} ));

	t_vec Pin = tl2::quat_vec_prod(quatB, _Pin);
	t_quat quatPin = tl2::rotation_quat(Pin, tl2::make_vec<t_vec>( {1, 0, 0} ));

	m_rotCoord = quatPin * quatB;
}


template<class t_real, class t_cplx, int ORDER>
void Skx<t_real, t_cplx, ORDER>::SetG(t_real h, t_real k, t_real l)
{
	m_Grlu = tl2::make_vec<t_vec>({h,k,l});

	t_vec G = m_Grlu / tl2::veclen(m_Grlu);
	G = tl2::quat_vec_prod(m_rotCoord, G);
	m_projNeutron = tl2::unit_m<t_mat>(3) - tl2::outer<t_vec, t_mat>(G, G);
}


/**
 * query the dispersion
 */
template<class t_real, class t_cplx, int ORDER>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
Skx<t_real, t_cplx, ORDER>::GetDisp(t_real h, t_real k, t_real l, t_real minE, t_real maxE) const
{
	t_vec Qrlu = tl2::make_vec<t_vec>( {h, k, l} );
	t_vec qrlu = Qrlu - m_Grlu;
	t_vec qkh = qrlu / G_KH_RLU_29K;

	qkh = tl2::quat_vec_prod(m_rotCoord, qkh);
	t_real _l = qkh[2];
	qkh.resize(2, true);

	t_vec Qmagrlu = tl2::prod_mv(m_Binv, qkh);
	t_vec Gmagrlu = tl2::make_vec<t_vec>({
		std::round(Qmagrlu[0]), std::round(Qmagrlu[1]) });


	static const std::vector<t_vec> sats =
	{
		tl2::make_vec<t_vec>({0, 0}),
		tl2::make_vec<t_vec>({0, -1}), tl2::make_vec<t_vec>({0, +1}),
		tl2::make_vec<t_vec>({-1, 0}), tl2::make_vec<t_vec>({+1, 0}),
		tl2::make_vec<t_vec>({-1, -1}), tl2::make_vec<t_vec>({+1, +1}),
	};

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

	return GetSpecWeights(int(Gmagrlu[0]), int(Gmagrlu[1]),
		qmagrlu[0], qmagrlu[1], _l, minE, maxE);
}
