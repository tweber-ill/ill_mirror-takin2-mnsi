/**
 * Free energy in the helimagnetic phase
 * @author tweber@ill.fr
 * @date mid-16, jul-18
 * @desc This file implements the theoretical helimagnon model by M. Garst and J. Waizner, references:
 *	- https://doi.org/10.1088/1361-6463/aa7573
 *	- https://kups.ub.uni-koeln.de/7937/
 *	- https://doi.org/10.1103/PhysRevLett.115.097203
 *	- personal communications with M. Garst, 2014-2019.
 * @desc This file is based on:
 *	- the descriptions and Mathematica implementations of the different helimagnon model versions by M. Garst and J. Waizner, 2014-2018,
 *	- the 2015 and 2016 Python implementations by G. Brandl and M. Kugler of the first version of the helimagnon model.
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
	// #1
	std::vector<int> m_idx1;
	m_idx1.resize(ORDER);
	std::iota(m_idx1.begin(), m_idx1.end(), 1);

	// #2
	for(int i=0; i<ORDER*2+1; ++i)
		std::copy(m_idx1.begin(), m_idx1.end(), std::back_inserter(m_idx2[0]));

	for(int i=-ORDER; i<=ORDER; ++i)
		for(int j=0; j<ORDER; ++j)
			m_idx2[1].push_back(i);

	// #3
	for(int i=0; i<ORDER*2+1; ++i)
	{
		std::copy(m_idx2[0].begin(), m_idx2[0].end(), std::back_inserter(m_idx3[0]));
		std::copy(m_idx2[1].begin(), m_idx2[1].end(), std::back_inserter(m_idx3[1]));
	}

	for(int i=-ORDER; i<=ORDER; ++i)
		for(int j=0; j<(2*ORDER+1)*ORDER; ++j)
			m_idx3[2].push_back(i);

	// reduce #2
	decltype(m_idx2) idx2;
	for(std::size_t i=0; i<m_idx2[0].size(); ++i)
	{
		int val = -m_idx2[0][i] - m_idx2[1][i];
		if(std::abs(val) > ORDER)
		{ continue; }
		else
		{
			idx2[0].push_back(m_idx2[0][i]);
			idx2[1].push_back(m_idx2[1][i]);
			idx2[2].push_back(val);
		}
	}
	for(int i=0; i<3; ++i)
		m_idx2[i] = std::move(idx2[i]);

	// reduce #3
	decltype(m_idx3) idx3;
	for(std::size_t i=0; i<m_idx3[0].size(); ++i)
	{
		int val = -m_idx3[0][i] - m_idx3[1][i] - m_idx3[2][i];
		if(std::abs(val) > ORDER)
		{ continue; }
		else
		{
			idx3[0].push_back(m_idx3[0][i]);
			idx3[1].push_back(m_idx3[1][i]);
			idx3[2].push_back(m_idx3[2][i]);
			idx3[3].push_back(val);
		}
	}
	for(int i=0; i<4; ++i)
		m_idx3[i] = std::move(idx3[i]);
}


/**
 * free energy
 * e.g.: M = (1+i, 1-i, 0), prop_1 = (0, 0, kh), prop_2 = (0, 0, -i*kh)
 *
 * Heisenberg term:
 *	- discrete, real: -J_ij S_i * S_j
 *	- continuous, real: J (del M(r))^2, Fourier: J k^2 |M(k)|^2
 *
 * DM term:
 *	- discrete, real: D_ij * (S_i cross S_j)
 *	- continuous, real: D M(r) (del cross M(r))
 */
template<class t_real, class t_cplx, int ORDER>
t_real Heli<t_real, t_cplx, ORDER>::F()
{
	// only z component used
	auto fourier0 = m_fourier[0][2];

	auto fourier_no0 = m_fourier;
	fourier_no0.erase(fourier_no0.begin());

	// add complex conjugate
	auto fourier_full = m_fourier;
	for(std::size_t i=m_fourier.size()-1; i>=1; --i)
		fourier_full.push_back(tl2::conjugate_vec(m_fourier[i]));


	std::vector<t_cplx> fourier_full_comp[3];
	for(const auto &vecC : fourier_full)
	{
		auto tup = split_vec3d(vecC);
		fourier_full_comp[0].push_back(std::get<0>(tup));
		fourier_full_comp[1].push_back(std::get<1>(tup));
		fourier_full_comp[2].push_back(std::get<2>(tup));
	}

	// lengths
	std::vector<t_real> lens;
	for(std::size_t i=0; i<fourier_no0.size(); ++i)
	{
		const auto [fourier_real, fourier_imag] = tl2::split_cplx_vec<
			ublas::vector<std::complex<t_real>>, ublas::vector<t_real>>(fourier_no0[i]);

		t_real len = ublas::inner_prod(fourier_real, fourier_real) +
			ublas::inner_prod(fourier_imag, fourier_imag);
		lens.push_back(2. * len);
	}


	// free energy
	t_cplx cF = 0;

	// dip
	cF += g_chi<t_real>/3. * fourier0*fourier0;
	for(std::size_t i=0; i<fourier_no0.size(); ++i)
		cF += 2. * g_chi<t_real> * std::norm(fourier_no0[i][2]);

	// dmi
	for(std::size_t i=0; i<ORDER; ++i)
	{
		t_real k = t_real(i+1);

		cF += 8. * m_pitch * k *
			( fourier_no0[i][0].real() * fourier_no0[i][1].imag() -
			fourier_no0[i][1].real() * fourier_no0[i][0].imag() );
	}

	// hoc
	for(std::size_t i=0; i<ORDER; ++i)
	{
		t_real k = t_real(i+1);

		cF += g_hoc<t_real> * m_pitch*m_pitch * k*k*k*k * lens[i];
	}

	// phi^2 & phi^4
	cF += (m_T + 1.) * fourier0*fourier0  +  fourier0*fourier0*fourier0*fourier0;
	for(std::size_t i=0; i<ORDER; ++i)
	{
		t_real k = t_real(i+1);

		// Heisenberg
		cF += m_pitch*m_pitch * k*k * lens[i];

		cF += (m_T + 1. + fourier0*fourier0) * lens[i];
	}

	for(std::size_t i=0; i<m_idx2[0].size(); ++i)
	{
		cF += 2. * fourier0 * get_comp(fourier_full_comp[2], m_idx2[0][i]) *
		(
			get_comp(fourier_full_comp[0], m_idx2[1][i]) * get_comp(fourier_full_comp[0], m_idx2[2][i]) +
			get_comp(fourier_full_comp[1], m_idx2[1][i]) * get_comp(fourier_full_comp[1], m_idx2[2][i]) +
			get_comp(fourier_full_comp[2], m_idx2[1][i]) * get_comp(fourier_full_comp[2], m_idx2[2][i])
		);
	}

	// phi^4
	for(std::size_t i=0; i<m_idx3[0].size(); ++i)
	{
		cF += 2. *
		(
			get_comp(fourier_full_comp[0], m_idx3[0][i]) * get_comp(fourier_full_comp[0], m_idx3[1][i]) +
			get_comp(fourier_full_comp[1], m_idx3[0][i]) * get_comp(fourier_full_comp[1], m_idx3[1][i]) +
			get_comp(fourier_full_comp[2], m_idx3[0][i]) * get_comp(fourier_full_comp[2], m_idx3[1][i])
		) *
		(
			get_comp(fourier_full_comp[0], m_idx3[2][i]) * get_comp(fourier_full_comp[0], m_idx3[3][i]) +
			get_comp(fourier_full_comp[1], m_idx3[2][i]) * get_comp(fourier_full_comp[1], m_idx3[3][i]) +
			get_comp(fourier_full_comp[2], m_idx3[2][i]) * get_comp(fourier_full_comp[2], m_idx3[3][i])
		);
	}


	// zee
	cF += -m_B*fourier0;
	return cF.real();
}


/**
 * set fourier components
 */
template<class t_real, class t_cplx, int ORDER>
void Heli<t_real, t_cplx, ORDER>::SetFourier(const std::vector<ublas::vector<t_cplx>> &fourier)
{
	m_fourier = fourier;
	while(m_fourier.size() < ORDER_FOURIER+1)
		m_fourier.emplace_back(tl2::make_vec<ublas::vector<t_cplx>>({0., 0., 0.}));
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
	const t_real sqrtfac = std::sqrt(Brel*Brel + 2.*Brel2*Brel2);
	const t_real E_scale_fac = (Brel*Brel + 2.*Brel2*Brel2) / std::sqrt(Brel2*Brel2) * g_g<t_real>*g_muB<t_real>*m_Bc2 * Brel2;


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
		sqrtfac, E_scale_fac,
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
