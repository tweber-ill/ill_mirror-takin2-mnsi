/**
 * MnSi constants
 * @author Tobias Weber <tweber@ill.fr>
 * @date jul-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#ifndef __MNSI_CONSTS_H__
#define __MNSI_CONSTS_H__

#include "tlibs2/libs/math17.h"
#include "tlibs2/libs/units.h"


/**
 * constants
 * The values are from the theoretical models by. M. Garst and J. Waizner:
 *	- https://doi.org/10.1088/1361-6463/aa7573
 *	- https://kups.ub.uni-koeln.de/7937/
 *	- Personal communications with M. Garst, 2017-2020.
 * Further values, e.g. g_hoc and g_chi, are from the paper and its Python code:
 *	- https://doi.org/10.1103/PhysRevLett.115.097203
 */
template<class t_real = double> constexpr t_real g_pi = tl2::pi<t_real>;
template<class t_real = double> constexpr t_real g_kB
	= static_cast<t_real>(tl2::kB<t_real>/tl2::meV<t_real>*tl2::kelvin<t_real>);
template<class t_real = double> constexpr t_real g_muB
	= static_cast<t_real>(tl2::muB<t_real>/tl2::meV<t_real>*tl2::tesla<t_real>);
template<class t_real = double> constexpr t_real g_chi = 0.34;
template<class t_real = double> constexpr t_real g_hoc = -0.0073;
template<class t_real = double> constexpr t_real g_a = 4.558;
template<class t_real = double> constexpr t_real g_g = 2.;
template<class t_real = double> constexpr t_real g_kh_A_29K = 0.039;
template<class t_real = double> constexpr t_real g_kh_rlu_29K = g_kh_A_29K<t_real> / (t_real(2.)*g_pi<t_real> / g_a<t_real>);


/**
 * helix pitch in inverse angstroms
 */
template<class t_real = double>
constexpr t_real g_kh_A(t_real T)
{
	constexpr t_real T1 = 20.;
	constexpr t_real T2 = 29.;

	constexpr t_real A1 = 0.036;			// at 20 K
	constexpr t_real A2 = g_kh_A_29K<t_real>;	// at 29 K

	constexpr t_real m = (A2-A1) / (T2-T1);
	constexpr t_real t = A2 - m*T2;

	return m*T + t;
}


/**
 * helix pitch in rlu
 */
template<class t_real = double>
constexpr t_real g_kh_rlu(t_real T)
{
	return g_kh_A<t_real>(T) / (t_real(2.)*g_pi<t_real> / g_a<t_real>);
}


/**
 * get upper critical field
 *	- Theoretical values calculated with "heli.dat_Bc2.gpl" script generated by heliphase.cpp
 *	- Experimental values were measured and provided by A. Bauer, 2015.
 */
template<class t_real=double>
t_real get_bc2(t_real T, bool use_theo_units=1)
{
	if(use_theo_units)
	{
		//const t_real amp = 1.62567;		// scaling
		//const t_real ex = 0.46491;		// critical exponent
		const t_real amp = 1.42925;		// scaling
		const t_real ex = 0.499501;		// critical exponent

		if(T >= 0.) return 0.;
		return amp * std::pow(-T, ex);

		//return g_g<t_real> * std::sqrt(-0.5 - 0.5*T);
	}
	else
	{
		t_real p1[] = { 5.63617388e-01, 5.52570801e-02,  2.20736277e+01, 7.39287474e-02, -5.32767610e-04 };
		t_real p2[] = { 0.07075453, -0.08217821, 30.00000534, 9.19469462,  0.38951838 };

		const t_real *p = (T<=20 ? p1 : p2);
		t_real Tc = p[2];
		t_real tau = (Tc-T) / Tc;

		if(T >= Tc) return 0.;
		return p[0]*std::pow(tau, p[1]) * (1. + p[3]*std::pow(tau, p[4]));
	}
}


#endif
