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
 *	- The 2016 Python implementations by M. Kugler and G. Brandl of the first version of the skyrmion model.
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
	// B matrix
	tl2::inverse(m_Bmat, m_Binv);

	// rotation in lab
	t_mat curRot = tl2::unit_m<t_mat>(2);
	const auto rot60 = tl2::rotation_matrix_2d<t_mat>(tl2::d2r<t_real>(60));

	//m_rot.push_back(curRot);
	for(int i=0; i<6; ++i)
	{
		curRot = tl2::prod_mm(curRot, rot60);
		m_rot.push_back(curRot);
	}

	// rotation in rlu
	for(const auto& rot : m_rot)
	{
		auto M = tl2::transform<t_mat>(rot, m_Bmat, 0);
		m_rotRlu.emplace_back(M);
	}


	// all peaks
	for(t_real h=-ORDER; h<ORDER+1; ++h)
		for(t_real k=-ORDER; k<ORDER+1; ++k)
			if(std::abs(h-k) <= t_real(ORDER))
				m_allpeaks.emplace_back(tl2::make_vec<t_vec>({h,k}));

	for(t_real h=0; h<ORDER+1; ++h)
		for(t_real k=0; k<ORDER; ++k)
			if(h>k)
			{
				t_vec vec_rlu = tl2::make_vec<t_vec>({h,k});
				t_vec vec_lab = tl2::prod_mv(m_Bmat, vec_rlu);
				vec_lab /= tl2::veclen(vec_lab);

				m_peaks_60.emplace_back(std::move(vec_rlu));
				m_peaks_60_lab.emplace_back(std::move(vec_lab));
			}


	// the top and bottom 180 degrees
	std::vector<t_vec> peaks_180_t, peaks_180_b;

	for(const auto& rot : m_rotRlu)
	{
		std::vector<t_vec> _peaks;
		for(const auto &vec : m_peaks_60)
		{
			t_vec pk_rlu = tl2::prod_mv(rot, vec);
			t_vec pk = tl2::prod_mv(m_Bmat, pk_rlu);

			if(tl2::float_equal<t_real>(pk[1], 0.))
			{
				if(pk[0] < 0.)
					peaks_180_b.push_back(pk_rlu);
				else
					peaks_180_t.push_back(pk_rlu);
			}
			else
			{
				if(pk[1] < 0.)
					peaks_180_b.push_back(pk_rlu);
				else
					peaks_180_t.push_back(pk_rlu);
			}

			_peaks.emplace_back(pk_rlu);
		}
		m_peaks_360.emplace_back(_peaks);
	}


	// #1
	for(const auto& vec : peaks_180_t)
		m_idx1.emplace_back(std::make_pair(
			int(std::round(vec[0])), int(std::round(vec[1]))));

	// #2
	for(std::size_t i=0; i<m_allpeaks.size(); ++i)
	{
		std::copy(m_idx1.begin(), m_idx1.end(), std::back_inserter(m_idx2[0]));
		for(std::size_t j=0; j<peaks_180_t.size(); ++j)
			m_idx2[1].emplace_back(std::make_pair(
				int(std::round(m_allpeaks[i][0])),
				int(std::round(m_allpeaks[i][1]))));
	}

	// #3
	for(std::size_t i=0; i<m_allpeaks.size(); ++i)
	{
		std::copy(m_idx2[0].begin(), m_idx2[0].end(), std::back_inserter(m_idx3[0]));
		std::copy(m_idx2[1].begin(), m_idx2[1].end(), std::back_inserter(m_idx3[1]));

		for(std::size_t j=0; j<peaks_180_t.size(); ++j)
			for(std::size_t k=0; k<m_allpeaks.size(); ++k)
				m_idx3[2].emplace_back(std::make_pair(
					int(std::round(m_allpeaks[i][0])),
					int(std::round(m_allpeaks[i][1]))));
	}


	// reduce #2
	decltype(m_idx2) idx2;
	for(std::size_t i=0; i<m_idx2[0].size(); ++i)
	{
		auto val1 = -m_idx2[0][i].first - m_idx2[1][i].first;
		auto val2 = -m_idx2[0][i].second - m_idx2[1][i].second;
		if(std::abs(val1) > ORDER || std::abs(val2) > ORDER || std::abs(val1-val2) > ORDER)
		{ continue; }
		else
		{
			idx2[0].push_back(m_idx2[0][i]);
			idx2[1].push_back(m_idx2[1][i]);
			idx2[2].emplace_back(std::make_pair(val1, val2));
		}
	}
	for(int i=0; i<3; ++i)
		m_idx2[i] = std::move(idx2[i]);

	// reduce #3
	decltype(m_idx3) idx3;
	for(std::size_t i=0; i<m_idx3[0].size(); ++i)
	{
		auto val1 = -m_idx3[0][i].first - m_idx3[1][i].first - m_idx3[2][i].first;
		auto val2 = -m_idx3[0][i].second - m_idx3[1][i].second - m_idx3[2][i].second;
		if(std::abs(val1) > ORDER || std::abs(val2) > ORDER || std::abs(val1-val2) > ORDER)
		{ continue; }
		else
		{
			idx3[0].push_back(m_idx3[0][i]);
			idx3[1].push_back(m_idx3[1][i]);
			idx3[2].push_back(m_idx3[2][i]);
			idx3[3].emplace_back(std::make_pair(val1, val2));
		}
	}
	for(int i=0; i<4; ++i)
		m_idx3[i] = std::move(idx3[i]);
}


template<class t_real, class t_cplx, int ORDER>
ublas::matrix<typename Skx<t_real, t_cplx, ORDER>::t_vec_cplx>
Skx<t_real, t_cplx, ORDER>::GetFullFourier() const
{
	constexpr auto imag = t_cplx(0,1);

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
				({ -vecPk[1]*m_fourier[ihx+1][0],
				   +vecPk[0]*m_fourier[ihx+1][0] });
			fourier = tl2::prod_mv(rot, fourier);

			get_comp(M, idx1, idx2) = tl2::make_vec<t_vec_cplx>
				({
					imag*fourier[0],
					imag*fourier[1],
					m_fourier[ihx+1][2]
				});
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
	auto fourier0 = tl2::veclen(m_M(0,0));

	// dip
	cF += g_chi<t_real>/3. * fourier0*fourier0;


	for(const auto& pair : m_idx1)
	{
		t_vec pos_rlu = tl2::make_vec<t_vec>({t_real(pair.first), t_real(pair.second)});
		t_vec pos_lab = tl2::prod_mv(m_Bmat, pos_rlu);

		const auto& m = get_comp(m_M, pair.first, pair.second);
		const auto len = tl2::inner_cplx(m, m);

		// dmi
		cF += 8. * m_pitch * imag*(pos_lab[0]*m[1]*m[2] - pos_lab[1]*m[2]*m[0]);

		// hoc
		cF += 2. * g_hoc<t_real> * m_pitch*m_pitch * len
			* std::pow(pos_rlu[0]*pos_rlu[0] + pos_rlu[1]*pos_rlu[1] - pos_rlu[0]*pos_rlu[1], 2.);

		// phi^2 & phi^4
		cF += 2. * m_pitch*m_pitch * len
			* (pos_rlu[0]*pos_rlu[0] + pos_rlu[1]*pos_rlu[1] - pos_rlu[0]*pos_rlu[1]);
		cF += (m_T + 1. + fourier0*fourier0) * 2.*len;
	}


	// phi^2 & phi^4
	cF += (m_T + 1.) * fourier0*fourier0  +  fourier0*fourier0*fourier0*fourier0;


	for(std::size_t i=0; i<m_idx2[0].size(); ++i)
	{
		const auto& m1 = get_comp(m_M, m_idx2[0][i].first, m_idx2[0][i].second);
		const auto& m2 = get_comp(m_M, m_idx2[1][i].first, m_idx2[1][i].second);
		const auto& m3 = get_comp(m_M, m_idx2[2][i].first, m_idx2[2][i].second);

		// the x and y components are purely imaginary, z purely real
		cF += 2. * fourier0 * m1[2] * tl2::inner(m2, m3);
	}


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
	cF += -m_B*fourier0;

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
	// indices

	// 2
	std::vector<std::pair<int, int>> idx2[3];
	for(std::size_t i=0; i<m_allpeaks.size(); ++i)
	{
		for(std::size_t j=0; j<m_allpeaks.size(); ++j)
		{
			idx2[0].emplace_back(std::make_pair(
				int(std::round(m_allpeaks[j][0])),
				int(std::round(m_allpeaks[j][1]))));
			idx2[1].emplace_back(std::make_pair(
				int(std::round(m_allpeaks[i][0])),
				int(std::round(m_allpeaks[i][1]))));
		}
	}

	// 3
	std::vector<std::pair<int, int>> idx3[4];
	for(std::size_t i=0; i<m_allpeaks.size(); ++i)
	{
		std::copy(idx2[0].begin(), idx2[0].end(), std::back_inserter(idx3[0]));
		std::copy(idx2[1].begin(), idx2[1].end(), std::back_inserter(idx3[1]));

		for(std::size_t j=0; j<m_allpeaks.size(); ++j)
			for(std::size_t k=0; k<m_allpeaks.size(); ++k)
				idx3[2].emplace_back(std::make_pair(
					int(std::round(m_allpeaks[i][0])),
					int(std::round(m_allpeaks[i][1]))));
	}


	// reduce 2
	decltype(idx2) _idx2;
	for(std::size_t i=0; i<idx2[0].size(); ++i)
	{
		auto val1 = idx2[0][i].first - idx2[1][i].first;
		auto val2 = idx2[0][i].second - idx2[1][i].second;
		if(std::abs(val1) > ORDER || std::abs(val2) > ORDER || std::abs(val1-val2) > ORDER)
		{
			continue;
		}
		else
		{
			_idx2[0].push_back(idx2[0][i]);
			_idx2[1].push_back(idx2[1][i]);
			_idx2[2].emplace_back(std::make_pair(val1, val2));
		}
	}
	for(int i=0; i<3; ++i)
		idx2[i] = std::move(_idx2[i]);

	// reduce 3
	decltype(idx3) _idx3;
	for(std::size_t i=0; i<idx3[0].size(); ++i)
	{
		auto val1 = idx3[0][i].first - idx3[1][i].first - idx3[2][i].first;
		auto val2 = idx3[0][i].second - idx3[1][i].second - idx3[2][i].second;
		if(std::abs(val1) > ORDER || std::abs(val2) > ORDER || std::abs(val1-val2) > ORDER)
		{
			continue;
		}
		else
		{
			_idx3[0].push_back(idx3[0][i]);
			_idx3[1].push_back(idx3[1][i]);
			_idx3[2].push_back(idx3[2][i]);
			_idx3[3].emplace_back(std::make_pair(val1, val2));
		}
	}
	for(int i=0; i<4; ++i)
		idx3[i] = std::move(_idx3[i]);
	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	// Mx_ij = eps_ikj M_k
	auto Mx = std::make_unique<std::array<t_mat_cplx, SIZE*SIZE*SIZE*SIZE>>();

	for(std::size_t i=0; i<idx2[2].size(); ++i)
	{
		const auto& vecM = get_comp(m_M, idx2[2][i].first, idx2[2][i].second);
		auto skew = tl2::skew<t_mat_cplx>(vecM);

		get_comp(*Mx, SIZE, idx2[0][i].first, idx2[0][i].second, idx2[1][i].first, idx2[1][i].second) = skew;
	}
	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	// fluct. F_q1q2_ij = 0.5 * d^2F/(dM_-q_1_q dM_q2_j)
	auto Fluc = std::make_unique<std::array<t_mat_cplx, SIZE*SIZE*SIZE*SIZE>>();

	for(std::size_t i=0; i<idx3[3].size(); ++i)
	{
		const auto& vecM1 = get_comp(m_M, idx3[2][i].first, idx3[2][i].second);
		const auto& vecM2 = get_comp(m_M, idx3[3][i].first, idx3[3][i].second);

		t_mat_cplx mat(3,3);
		for(int _midx1=0; _midx1<3; ++_midx1)
			for(int _midx2=0; _midx2<3; ++_midx2)
				mat(_midx1,_midx2) = vecM1[_midx1]*vecM2[_midx2];
		mat *= 4.;

		// delta term
		for(int d=0; d<3; ++d)
			mat(d,d) += 2.*tl2::inner(vecM1, vecM2);

		auto& oldmat = get_comp(*Fluc, SIZE, idx3[0][i].first, idx3[0][i].second, idx3[1][i].first, idx3[1][i].second);
		if(!oldmat.size1()) oldmat = 2.*mat; else oldmat += 2.*mat;
	}


	for(const auto &pk_rlu : m_allpeaks)
	{
		t_vec pk_lab = tl2::prod_mv(m_Bmat, pk_rlu);

		t_vec pos = pk_lab + q_lab;
		pos.resize(3, true);
		pos[2] = ql;

		t_real dipole = g_chi<t_real> / tl2::inner(pos, pos);

		t_mat_cplx mat(3,3);
		mat(0,0) = dipole * pos[0]*pos[0] + (1. + m_T) + tl2::inner(pos, pos);
		mat(0,1) = dipole * pos[0]*pos[1] - 2.*imag*pos[2];
		mat(0,2) = dipole * pos[0]*pos[2] + 2.*imag*pos[1];
		mat(1,0) = dipole * pos[1]*pos[0] + 2.*imag*pos[2];
		mat(1,1) = dipole * pos[1]*pos[1] + (1. + m_T) + tl2::inner(pos, pos);
		mat(1,2) = dipole * pos[1]*pos[2] - 2.*imag*pos[0];
		mat(2,0) = dipole * pos[2]*pos[0] - 2.*imag*pos[1];
		mat(2,1) = dipole * pos[2]*pos[1] + 2.*imag*pos[0];
		mat(2,2) = dipole * pos[2]*pos[2] + (1. + m_T) + tl2::inner(pos, pos);

		int idx1 = int(std::round(pk_rlu[0]));
		int idx2 = int(std::round(pk_rlu[1]));
		auto& oldmat = get_comp(*Fluc, SIZE, idx1, idx2, idx1, idx2);
		if(!oldmat.size1()) oldmat = 2.*mat; else oldmat += 2.*mat;
	}
	// ------------------------------------------------------------------------


	auto mk_2dim = [VIRTSIZE, CRYSSIZE, iGh, iGk](const /*auto&*/ decltype(*Mx)& arr) -> t_mat_cplx
	{
		std::vector<int> pks1(VIRTSIZE);
		std::iota(pks1.begin(), pks1.begin()+CRYSSIZE+1, 0);
		std::iota(pks1.begin()+CRYSSIZE+1, pks1.end(), -CRYSSIZE);

		std::vector<std::pair<int, int>> pks2;
		pks2.reserve(VIRTSIZE*VIRTSIZE);
		for(int k : pks1)
			for(int h : pks1)
				if(std::abs(h-k) <= CRYSSIZE)
					pks2.emplace_back(std::make_pair(h, k));

		auto mat = tl2::zero_matrix<t_mat_cplx>(3*pks2.size(), 3*pks2.size());
		for(std::size_t idx1=0; idx1<pks2.size(); ++idx1)
		{
			const auto& hk1 = pks2[idx1];
			for(std::size_t idx2=0; idx2<pks2.size(); ++idx2)
			{
				const auto& hk2 = pks2[idx2];

				const auto& comp = get_virt_comp(arr, SIZE, VIRTSIZE, ORDER,
					hk1.first+iGh, hk1.second+iGk, hk2.first+iGh, hk2.second+iGk);

				if(comp.size1() && comp.size2())
				{
					for(std::size_t x1=0; x1<3; ++x1)
						for(std::size_t x2=0; x2<3; ++x2)
							mat(idx1*3+x1, idx2*3+x2) = comp(x1,x2);
				}
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
	{ qh += epshkl; qk += epshkl; ql += epshkl; }

	ql = -ql;

	t_real sqrtfac = std::sqrt(t_real(-0.5 - m_T*0.5));
	// energy scaling depends on Hc2_int (and the sample's demagnetisation factor)
	// Hc2_int changes rapidly vs. T here at the border to the paramegnetic phase
	const t_real E_scale_fac_heli = 0.0387;		// calculated with heli.cpp
	t_real E_scale_fac = E_scale_fac_heli / (sqrtfac/10.);

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

	return GetSpecWeights(int(Gmagrlu[0]), int(Gmagrlu[1]), qmagrlu[0], qmagrlu[1], _l, minE, maxE);
}
