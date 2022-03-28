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


// instantiation
#ifdef DEF_HELI_ORDER
	#pragma message("Heli Order: " __TL2_STRCONV(DEF_HELI_ORDER))
	template class Heli<double, std::complex<double>, DEF_HELI_ORDER>;

	#ifdef __HACK_FULL_INST__
		template Heli<double, std::complex<double>, DEF_HELI_ORDER>::Heli();
		template void Heli<double, std::complex<double>, DEF_HELI_ORDER>::SetG(double, double, double);
		template void Heli<double, std::complex<double>, DEF_HELI_ORDER>::SetCoords(double, double, double);
	#endif
#endif


/**
 * constructor
 */
template<class t_real, class t_cplx, int ORDER>
Heli<t_real, t_cplx, ORDER>::Heli()
{
	for(int i=0; i<2; ++i)
		m_idx2[i].reserve(ORDER * (ORDER*2+1));
	for(int i=0; i<3; ++i)
		m_idx3[i].reserve(ORDER * (ORDER*2+1) * (ORDER*2+1));

	// unrolled indices for two loops
	for(int j=-ORDER; j<=ORDER; ++j)
	{
		for(int i=1; i<ORDER+1; ++i)
		{
			int k = -i-j;
			if(std::abs(k) > ORDER)
				continue;

			m_idx2[0].push_back(i);
			m_idx2[1].push_back(j);
			m_idx2[2].push_back(k);
		}
	}

	// unrolled indices for three loops
	for(int k=-ORDER; k<=ORDER; ++k)
	{
		for(int j=-ORDER; j<=ORDER; ++j)
		{
			for(int i=1; i<ORDER+1; ++i)
			{
				int l = -i-j-k;
				if(std::abs(l) > ORDER)
					continue;

				m_idx3[0].push_back(i);
				m_idx3[1].push_back(j);
				m_idx3[2].push_back(k);
				m_idx3[3].push_back(l);
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
	const auto& m = m_fourier;

	const auto& m0 = m[0];
	const auto m0_sq = tl2::inner_cplx(m0, m0);

	// free energy
	t_cplx cF = 0;

	// dip
	cF += g_chi<t_real>/3. * m0_sq;

	for(std::size_t i=1; i<ORDER+1; ++i)
	{
		const t_real q = t_real(i);
		const t_real q_sq = q*q;
		const auto m_sq = tl2::inner_cplx(m[i], m[i]);

		// dip
		cF += 2. * g_chi<t_real> * std::norm(m[i][2]);

		// dmi
		cF += 8. * q * (
			m[i][0].real() * m[i][1].imag() -
			m[i][1].real() * m[i][0].imag() );

		// hoc
		cF += 2. * g_hoc<t_real> * m_sq * q_sq*q_sq;

		// phi^2 & phi^4
		cF += 2. * m_sq * q_sq;
		cF += 2. * (m_T + 1. + m0_sq) * m_sq;
	}

	// phi^2 & phi^4
	cF += (m_T + 1.) * m0_sq + m0_sq*m0_sq;

	for(std::size_t i=0; i<m_idx2[0].size(); ++i)
	{
		const auto& m1 = get_comp(m, m_idx2[0][i]);
		const auto& m2 = get_comp(m, m_idx2[1][i]);
		const auto& m3 = get_comp(m, m_idx2[2][i]);

		cF += 2. * tl2::inner(m0, m1) * tl2::inner(m2, m3);
	}

	// phi^4
	for(std::size_t i=0; i<m_idx3[0].size(); ++i)
	{
		const auto& m1 = get_comp(m, m_idx3[0][i]);
		const auto& m2 = get_comp(m, m_idx3[1][i]);
		const auto& m3 = get_comp(m, m_idx3[2][i]);
		const auto& m4 = get_comp(m, m_idx3[3][i]);

		cF += 2. * tl2::inner(m1, m2) * tl2::inner(m3, m4);
	}

	// zee
	cF += -m_B * std::sqrt(m0_sq);
	return cF.real();
}


/**
 * set fourier components
 */
template<class t_real, class t_cplx, int ORDER>
void Heli<t_real, t_cplx, ORDER>::SetFourier(const std::vector<t_vec_cplx> &fourier)
{
	m_fourier = fourier;

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
	constexpr int SIZE = 2*ORDER+1;
	constexpr t_cplx imag = t_cplx(0,1);
	static const std::vector<t_real> empty;
	const t_real epshkl = 1e-5;

	t_vec qvec = tl2::make_vec<t_vec>({ qh, qk, -ql });

	if(tl2::float_equal<t_real>(qvec[2], 0., epshkl))
		qvec[2] += epshkl;


	const t_real Brel = m_B / m_Bc2;
	const t_real Brel2 = std::sqrt(0.5 - 0.5*Brel*Brel);

	// static susceptibility
	t_mat_cplx fluct = tl2::zero_m<t_mat_cplx>(3*SIZE, 3*SIZE);

	constexpr t_real A = g_hoc<t_real>;
	constexpr t_real A2 = A*A;
	constexpr t_real A3 = A2*A;

	static const/*expr*/ t_real _c1 = (-2.*imag*std::pow(2., 2./3.) * std::pow(3., 5./6.) * A *
		std::pow(-9.*A2 + std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))), 1./3.) /
		(std::pow(2., 1./3.) * (3.+imag*std::sqrt(3.)) * A + std::pow(3., 1./6.) * (std::sqrt(3.)-imag) *
		std::pow(-9.*A2 + std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))), 2./3.))).real();

	static const/*expr*/ t_real _c2 = (std::pow(std::pow(2., 1./3.) * (std::sqrt(3.)-3.*imag) * A -
		imag*std::pow(3.,1./6.)*(std::sqrt(3.)-imag) *
		std::pow(-9.*A2 + std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))), 2./3.), 2.) /
		(24.*std::pow(2., 1./3.) * A*std::pow(-27.*A2+3.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))), 2./3.))).real();

	static const/*expr*/ t_real _c3 = ((324.*A3 - 9.*imag*(std::sqrt(3.)-imag)*A2*std::pow(-54.*A2+6.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))), 1./3.) +
		(3.*imag + std::sqrt(3)) * std::sqrt(t_cplx(A3*(2.+27.*A))) * std::pow(-54.*A2+6.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))), 1./3.) -
		A*(36.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))) -3.*imag*std::pow(3., 1./6.) * std::pow(-18.*A2 + 2.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))), 2./3.) +
		std::pow(-54.*A2 + 6.*std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))), 2./3.))) /
		(2.*std::pow(2., 1./3.) * std::pow(3., 1./6.) * std::pow(-9.*A2+std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))), 2./3.) *
		(std::pow(2., 1./3.) * (-3.*imag+std::sqrt(3.))*A - imag*std::pow(3., 1./6) *
		(-imag+std::sqrt(3.)) * std::pow(-9.*A2+std::sqrt(t_cplx(3.*(A3*(2.+27.*A)))), 2./3.)))).real();

	const t_real _c4 = _c3 + 88./3. * Brel2*Brel2;
	const t_real _c5 = _c3 + 88./3. * Brel*Brel;

	for(int pos=0; pos<int(SIZE); ++pos)
	{
		t_real qz = qvec[2] - t_real(ORDER) + t_real(pos);
		const t_real q2 = qvec[0]*qvec[0] + qvec[1]*qvec[1] + qz*qz;

		fluct(3*pos + 0, 3*pos + 0) = (g_chi<t_real>*(qvec[0]*qvec[0] + qvec[1]*qvec[1]) + 2.*q2*(2.*_c1*qz + q2 + _c2*q2*q2 + _c4)) / (2.*q2);
		fluct(3*pos + 0, 3*pos + 1) = (g_chi<t_real>*(qvec[0] - imag*qvec[1])*(qvec[0] - imag*qvec[1])) / (2.*q2);
		fluct(3*pos + 0, 3*pos + 2) = ((qvec[0] - imag*qvec[1]) * (g_chi<t_real>*qz - 2*_c1*q2) * std::sqrt(2)) / (2.*q2);

		fluct(3*pos + 1, 3*pos + 0) = std::conj(fluct(3*pos + 0, 3*pos + 1));
		fluct(3*pos + 1, 3*pos + 1) = (g_chi<t_real>*(qvec[0]*qvec[0] + qvec[1]*qvec[1]) + 2.*q2*(-2.*_c1*qz + q2 + _c2*q2*q2 + _c4)) / (2.*q2);
		fluct(3*pos + 1, 3*pos + 2) = ((qvec[0] + imag*qvec[1]) * (g_chi<t_real>*qz + 2.*_c1*q2) * std::sqrt(2.)) / (2.*q2);

		fluct(3*pos + 2, 3*pos + 0) = std::conj(fluct(3*pos + 0, 3*pos + 2));
		fluct(3*pos + 2, 3*pos + 1) = std::conj(fluct(3*pos + 1, 3*pos + 2));
		fluct(3*pos + 2, 3*pos + 2) = (2.*g_chi<t_real>*qz*qz + 2.*q2*(q2 + _c2*q2*q2 + _c5)) / (2.*q2);

		if(pos > 1)
			fluct(3*pos + 1, 3*(pos-2) + 0) = Brel2 * Brel2 * 88./3.;
		if(pos > 0)
			fluct(3*pos + 1, 3*(pos-1) + 2) = fluct(3*pos + 2, 3*(pos-1) + 0) = Brel * Brel2 * 88./3.;
		if(pos < int(SIZE)-1)
			fluct(3*pos + 0, 3*(pos+1) + 2) = fluct(3*pos + 2, 3*(pos+1) + 1) = Brel * Brel2 * 88./3.;
		if(pos < int(SIZE)-2)
			fluct(3*pos + 0, 3*(pos+2) + 1) = Brel2 * Brel2 * 88./3.;
	}


	// Mx
	t_mat_cplx Mx = tl2::zero_m<t_mat_cplx>(3*SIZE, 3*SIZE);
	for(int pos=0; pos<int(SIZE); ++pos)
	{
		if(pos > 0)
		{
			Mx(3*pos + 2, 3*(pos-1) + 0) = imag*Brel2;
			Mx(3*pos + 1, 3*(pos-1) + 2) = -imag*Brel2;
		}

		Mx(3*pos + 0, 3*pos + 0) = -imag*Brel;
		Mx(3*pos + 1, 3*pos + 1) = imag*Brel;

		if(pos < int(SIZE)-1)
		{
			Mx(3*pos + 0, 3*(pos+1) + 2) = imag*Brel2;
			Mx(3*pos + 2, 3*(pos+1) + 1) = -imag*Brel2;
		}
	}


	// energies and weights
	return calc_weights<t_mat_cplx, t_vec_cplx, t_cplx, t_real>(
		Mx, fluct,
		m_bProjNeutron, m_projNeutron, m_polMat,
		1., g_g<t_real>*g_muB<t_real>*m_Bc2, // E scale factor
		minE, maxE,
		m_eveps, /*m_evlimit*/ -1., m_weighteps,
		m_filterzeroweight, m_onlymode,
		3*ORDER);
}


/**
 * set lattice vector and orthogonal projector
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
	{
		m_projNeutron = tl2::conjugate_mat(m_projNeutron);
	}
}


/**
 * rotates field to internal [001] convention
 */
template<class t_real, class t_cplx, int ORDER>
void Heli<t_real, t_cplx, ORDER>::SetCoords(t_real Bx, t_real By, t_real Bz)
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
