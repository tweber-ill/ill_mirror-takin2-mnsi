/**
 * Longitudinal fluctuations
 * @author Tobias Weber <tweber@ill.fr>
 * @date 13-jul-20
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "longfluct.h"
#include "constants.h"


Longfluct::Longfluct()
{
	SetPinning(1,-1,0, 0,0,1);
}


/**
 * generate skx satellites depending on the pinning and the up (i.e. B) vector
 */
void Longfluct::SetPinning(t_real Px, t_real Py, t_real Pz,
	t_real upx, t_real upy, t_real upz)
{
	m_up = tl2::make_vec<t_vec>({upx, upy, upz});
	m_up /= tl2::veclen(m_up);

	t_vec sat = tl2::make_vec<t_vec>({Px, Py, Pz});
	sat = sat / tl2::veclen(sat) * g_kh_rlu_29K<t_real>;

	t_mat rot = tl2::rotation_matrix(m_up, tl2::d2r(60.));

	m_sats.clear();
	m_sats.push_back(sat);
	for(int i=1; i<6; ++i)
		m_sats.emplace_back(tl2::prod_mv(rot, m_sats[i-1]));
}


/**
 * calculates the dynamical structure factor for the longitudinal fluctuations
 * using the formula from M. Garst, 10/jul/2020.
 */
Longfluct::t_real Longfluct::S_para(const Longfluct::t_vec& q, Longfluct::t_real E) const
{
	if(tl2::float_equal<t_real>(m_A, 0))
		return 0.;

	constexpr t_real Ec2 = 0.04;
	constexpr t_real kh = g_kh_rlu_29K<t_real>;

	t_real S = 0.;
	for(const t_vec& sat : m_sats)
	{
		S += g_kB<t_real> * m_T * m_A /
		(
			(m_inv_correl*m_inv_correl + tl2::inner(q-sat, q-sat)) / (kh*kh) +
			std::pow(m_Gamma*E / Ec2, 2.)
		);
	}
	return S;
}
