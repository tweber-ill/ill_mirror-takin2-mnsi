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


// constants
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


template<class t_real = double>
constexpr t_real g_kh_rlu(t_real T)
{
	return g_kh_A<t_real>(T) / (t_real(2.)*g_pi<t_real> / g_a<t_real>);
}


#define G_CHI g_chi<t_real>
#define G_HOC g_hoc<t_real>
#define G_G g_g<t_real>
#define G_MUB g_muB<t_real>
#define G_KH_RLU_29K g_kh_rlu_29K<t_real>


#endif
