/**
 * Longitudinal fluctuations
 * @author Tobias Weber <tweber@ill.fr>
 * @date 13-jul-20
 * @license GPLv2 (see 'LICENSE' file)
 */

#define USE_LAPACK
//#include "tlibs2/libs/math17.h"
#include "tlibs2-extras/math17.h"

#include <vector>


class Longfluct
{
public:
	using t_real = double;
	using t_vec = ublas::vector<t_real>;
	using t_mat = ublas::matrix<t_real>;

public:
	Longfluct();

	void SetPinning(t_real Px=1, t_real Py=-1, t_real Pz=0,
		t_real upx=0, t_real upy=0, t_real upz=1);

	t_real S_para(const t_vec& q, t_real E) const;

	void SetT(t_real T) { m_T = T; }
	void SetGamma(t_real G) { m_Gamma = G; }
	void SetInvCorrel(t_real k) { m_inv_correl = k; }
	void SetA(t_real A) { m_A = A; }

	t_real GetT() const { return m_T; }
	t_real GetGamma() const { return m_Gamma; }
	t_real GetInvCorrel() const { return m_inv_correl; }
	t_real GetA() const { return m_A; }

private:
	// skx satellites
	std::vector<t_vec> m_sats;

	t_vec m_up;

	t_real m_T = 29.;
	t_real m_Gamma = 0.5;
	t_real m_inv_correl = 0.1;
	t_real m_A = 0.;
};
