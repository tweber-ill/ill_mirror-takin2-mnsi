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
 *	- The 2015 and 2016 Python optimised implementations by G. Brandl and M. Kugler of the first version of the helimagnon model.
 *	  This present version started as a C++ port of that Python implementation by G. Brandl and M. Kugler,
 *	  that was then adapted to new theoretical model revisions provided by M. Garst.
 *	- The 2016 optimised Python implementations by M. Kugler and G. Brandl of the first version of the skyrmion model,
 *	  that also included improvements to the helimagnon code, which we ported to C++ in the present version.
 * @license GPLv2 (see 'LICENSE' file)
 */

#include <cmath>
#include <numeric>
#include <functional>
#include <iostream>

#include "heli.h"
#include "heli_inst.cxx"


/**
 * constructor
 */
template<class t_real, class t_cplx, int ORDER>
Heli<t_real, t_cplx, ORDER>::Heli()
{
	for(int i=0; i<2; ++i)
		m_idx2[i].reserve(SIZE * ORDER);
	for(int i=0; i<3; ++i)
		m_idx3[i].reserve(SIZE * SIZE * ORDER);

	for(int k=-ORDER; k<=ORDER; ++k)
	{
		for(int j=-ORDER; j<=ORDER; ++j)
		{
			for(int i=1; i<=ORDER; ++i)
			{
				// unrolled indices for three loops
				int l = - i - j - k;
				if(std::abs(l) <= ORDER)
				{
					m_idx3[0].push_back(i);
					m_idx3[1].push_back(j);
					m_idx3[2].push_back(k);
					m_idx3[3].push_back(l);
				}

				// unrolled indices for two loops
				int m = - i - j;
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
 * free energy
 */
template<class t_real, class t_cplx, int ORDER>
t_real Heli<t_real, t_cplx, ORDER>::F()
{
	const t_vec_cplx& m0 = m_fourier[0];
	const auto m0_sq = tl2::inner(m0, m0);

	// dipolar interaction
	const t_mat_cplx demag = tl2::diag_matrix<t_mat_cplx>({1./3., 1./3., 1./3.});
	t_cplx cF = g_chi<t_real> * tl2::inner(m0, tl2::prod_mv(demag, m0));

	// phi^2
	cF += (m_T_theo + 1.) * m0_sq;

	// phi^4
	cF += m0_sq*m0_sq;

	for(std::size_t i=1; i<=ORDER; ++i)
	{
		const t_vec_cplx& mi = m_fourier[i];
		t_vec_cplx mj = tl2::conjugate_vec(mi);
		const auto m_sq = tl2::inner(mi, mj);

		const t_real q = t_real(i); // q_vec = [0, 0, q]
		const t_real q_sq = q*q;

		// dipolar interaction
		cF += 2. * g_chi<t_real> * mi[2]*mj[2];

		// dmi
		cF += -4. * t_cplx(0, 1) * (-mi[0]*q*mj[1] + mi[1]*q*mj[0]);

		// phi^2
		cF += 2. * m_sq * q_sq;
		cF += 2. * (m_T_theo + 1.) * m_sq;

		// phi^4
		cF += 2. * m0_sq * m_sq;

		// high-order correction
		cF += 2. * g_hoc<t_real> * m_sq * q_sq*q_sq;
	}

	// phi^4
	for(std::size_t i=0; i<m_idx2[0].size(); ++i)
	{
		const auto& m1 = get_comp(m_fourier, m_idx2[0][i]);
		const auto& m2 = get_comp(m_fourier, m_idx2[1][i]);
		const auto& m3 = get_comp(m_fourier, m_idx2[2][i]);

		cF += 2. * tl2::inner(m0, m1) * tl2::inner(m2, m3);
	}
	for(std::size_t i=0; i<m_idx3[0].size(); ++i)
	{
		const auto& m1 = get_comp(m_fourier, m_idx3[0][i]);
		const auto& m2 = get_comp(m_fourier, m_idx3[1][i]);
		const auto& m3 = get_comp(m_fourier, m_idx3[2][i]);
		const auto& m4 = get_comp(m_fourier, m_idx3[3][i]);

		cF += 2. * tl2::inner(m1, m2) * tl2::inner(m3, m4);
	}

	// zeeman shift
	cF += -m_B_theo * std::sqrt(m0_sq);
	return cF.real();
}


/**
 * set fourier components
 */
template<class t_real, class t_cplx, int ORDER>
void Heli<t_real, t_cplx, ORDER>::SetFourier(const std::vector<t_vec_cplx> &fourier, bool symm)
{
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

	// add complex conjugate bottom peaks
	for(std::size_t i=ORDER_FOURIER; i>=1; --i)
		m_fourier.push_back(tl2::conjugate_vec(m_fourier[i]));
}


/**
 * energies and spectral weights
 */
template<class t_real, class t_cplx, int ORDER>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
Heli<t_real, t_cplx, ORDER>::GetSpecWeights(t_real qh, t_real qk, t_real ql, t_real minE, t_real maxE) const
{
	constexpr t_cplx imag = t_cplx(0, 1);

	avoid_G(qh, qk, ql, m_eps);
	const t_vec qvec = tl2::make_vec<t_vec>({ qh, qk, -ql });
	const t_cplx hx = qvec[0] - imag*qvec[1];

	constexpr t_real A1 = g_hoc<t_real>;
	constexpr t_real A2 = A1*A1;
	constexpr t_real A3 = A2*A1;

	static const t_real _c1 = (-2.*imag*std::pow(2., 2./3.) * std::pow(3., 5./6.) * A1 *
		std::pow(-9.*A2 + std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 1./3.) /
		(std::pow(2., 1./3.) * (3.+imag*std::sqrt(3.)) * A1 + std::pow(3., 1./6.) * (std::sqrt(3.)-imag) *
		std::pow(-9.*A2 + std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.))).real();
	static const t_real _c2 = (std::pow(std::pow(2., 1./3.) * (std::sqrt(3.)-3.*imag) * A1 -
		imag*std::pow(3.,1./6.)*(std::sqrt(3.)-imag) *
		std::pow(-9.*A2 + std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.), 2.) /
		(24.*std::pow(2., 1./3.) * A1*std::pow(-27.*A2+3.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.))).real();
	static const t_real _c3 = ((324.*A3 - 9.*imag*(std::sqrt(3.)-imag)*A2*std::pow(-54.*A2+6.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 1./3.) +
		(3.*imag + std::sqrt(3)) * std::sqrt(t_cplx(A3*(2.+27.*A1))) * std::pow(-54.*A2+6.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 1./3.) -
		A1*(36.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))) -3.*imag*std::pow(3., 1./6.) * std::pow(-18.*A2 + 2.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.) +
		std::pow(-54.*A2 + 6.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.))) /
		(2.*std::pow(2., 1./3.) * std::pow(3., 1./6.) * std::pow(-9.*A2+std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.) *
		(std::pow(2., 1./3.) * (-3.*imag+std::sqrt(3.))*A1 - imag*std::pow(3., 1./6) *
		(-imag+std::sqrt(3.)) * std::pow(-9.*A2+std::sqrt(t_cplx(3.*(A3*(2.+27.*A1)))), 2./3.)))).real();

	const t_real Brel = m_B / m_Bc2;
	const t_real Brel2 = std::sqrt(0.5 - 0.5*Brel*Brel);

	t_mat_cplx Mx = tl2::zero_m<t_mat_cplx>(3*SIZE, 3*SIZE);
	t_mat_cplx fluct = tl2::zero_m<t_mat_cplx>(3*SIZE, 3*SIZE);

	for(int pk=0; pk<SIZE; ++pk)
	{
		const t_real qz = qvec[2] + t_real(pk) - t_real(ORDER);
		const t_real qxy2 = qvec[0]*qvec[0] + qvec[1]*qvec[1];
		const t_real q2 = qxy2 + qz*qz;

		// M-cross matrix
		Mx(3*pk + 0, 3*pk + 0) = -imag*Brel;
		Mx(3*pk + 1, 3*pk + 1) = imag*Brel;

		if(pk > 0)
		{
			Mx(3*pk + 2, 3*(pk-1) + 0) = imag*Brel2;
			Mx(3*pk + 1, 3*(pk-1) + 2) = -imag*Brel2;
		}

		if(pk < SIZE-1)
		{
			Mx(3*pk + 0, 3*(pk+1) + 2) = imag*Brel2;
			Mx(3*pk + 2, 3*(pk+1) + 1) = -imag*Brel2;
		}

		// fluctuation matrix
		auto get_dip = [this](t_real q, t_real q_sq) -> t_real
		{
			if(tl2::float_equal<t_real>(q_sq, 0., m_eps))
				return 0.;
			return g_chi<t_real>/q_sq * q;
		};

		fluct(3*pk + 0, 3*pk + 0) = 0.5*get_dip(qxy2, q2) + 2.*_c1*qz + q2 + _c2*q2*q2 + _c3 + 88./3.*Brel2*Brel2;
		fluct(3*pk + 0, 3*pk + 1) = 0.5*get_dip(1., q2) * hx * hx;
		fluct(3*pk + 0, 3*pk + 2) = (0.5*get_dip(qz, q2) - _c1) * std::sqrt(2) * hx;

		fluct(3*pk + 1, 3*pk + 0) = std::conj(fluct(3*pk + 0, 3*pk + 1));
		fluct(3*pk + 1, 3*pk + 1) = fluct(3*pk + 0, 3*pk + 0) - 4.*_c1*qz;
		fluct(3*pk + 1, 3*pk + 2) = (0.5*get_dip(qz, q2) + _c1) * std::sqrt(2.) * std::conj(hx);

		fluct(3*pk + 2, 3*pk + 0) = std::conj(fluct(3*pk + 0, 3*pk + 2));
		fluct(3*pk + 2, 3*pk + 1) = std::conj(fluct(3*pk + 1, 3*pk + 2));
		fluct(3*pk + 2, 3*pk + 2) = get_dip(qz*qz, q2) + q2 + _c2*q2*q2 + _c3 + 88./3.*Brel*Brel;

		if(pk > 1)
			fluct(3*pk + 1, 3*(pk-2) + 0) = 88./3.*Brel2*Brel2;
		if(pk > 0)
			fluct(3*pk + 1, 3*(pk-1) + 2) = fluct(3*pk + 2, 3*(pk-1) + 0) = 88./3.*Brel*Brel2;
		if(pk < SIZE-1)
			fluct(3*pk + 0, 3*(pk+1) + 2) = fluct(3*pk + 2, 3*(pk+1) + 1) = 88./3.*Brel*Brel2;
		if(pk < SIZE-2)
			fluct(3*pk + 0, 3*(pk+2) + 1) = 88./3.*Brel2*Brel2;
	}

	// energies and weights
	return calc_weights<t_mat_cplx, t_vec_cplx, t_cplx, t_real>(
		Mx * g_g<t_real>, fluct,
		m_bProjNeutron, m_projNeutron, m_polMat,
		g_muB<t_real> * m_Bc2, g_g<t_real>,
		minE, maxE, m_eveps, /*m_evlimit*/ -1., m_weighteps,
		m_filterzeroweight, m_onlymode, 3*ORDER);
}


/**
 * set the lattice vector and orthogonal projector
 */
template<class t_real, class t_cplx, int ORDER>
void Heli<t_real, t_cplx, ORDER>::SetG(t_real h, t_real k, t_real l)
{
	bool bInChiralBase = true;
	m_Grlu = tl2::make_vec<t_vec>({h,k,l});

	t_vec _G = m_Grlu / tl2::veclen(m_Grlu);
	_G = tl2::quat_vec_prod(m_rotCoord, _G);
	t_vec_cplx G = _G;

	if(bInChiralBase)
	{
		t_mat_cplx chiral = get_chiralbasismat<t_mat_cplx, t_vec_cplx>();
		G = tl2::prod_mv<t_vec_cplx, t_mat_cplx>(chiral, G);
	}

	m_projNeutron = tl2::unit_m<t_mat_cplx>(3);
	m_projNeutron -= tl2::outer_cplx<t_vec_cplx, t_mat_cplx>(G, G);

	if(bInChiralBase)
		m_projNeutron = tl2::conjugate_mat(m_projNeutron);
}


/**
 * rotates field to internal [001] convention
 */
template<class t_real, class t_cplx, int ORDER>
void Heli<t_real, t_cplx, ORDER>::SetCoords(t_real Bx, t_real By, t_real Bz, t_real Px, t_real Py, t_real pZ)
{
	t_vec B = tl2::make_vec<t_vec>( {Bx, By, Bz} );
	t_quat quatB = tl2::rotation_quat(B, tl2::make_vec<t_vec>( {0, 0, 1} ));
	m_rotCoord = quatB;
}


/**
 * dispersion
 */
template<class t_real, class t_cplx, int ORDER>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
Heli<t_real, t_cplx, ORDER>::GetDisp(t_real h, t_real k, t_real l, t_real minE, t_real maxE) const
{
	t_vec Qvec = tl2::make_vec<t_vec>({ h, k, l });
	t_vec qvec = tl2::quat_vec_prod(m_rotCoord, Qvec) - tl2::quat_vec_prod(m_rotCoord, m_Grlu);
	qvec /= g_kh_rlu<t_real>(m_T);
	return GetSpecWeights(qvec[0], qvec[1], qvec[2], minE, maxE);
}
