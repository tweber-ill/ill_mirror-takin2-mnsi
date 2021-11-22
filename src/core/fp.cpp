/**
 * Field-polarised phase
 * @author Tobias Weber <tweber@ill.fr>
 * @date dec-16, sep-18
 * @desc This file implements the theoretical magnon model by M. Garst and J. Waizner, see:
 * 	- https://doi.org/10.1088/1361-6463/aa7573
 * 	- https://kups.ub.uni-koeln.de/7937/
 *	- personal communications with M. Garst, 2016-2019.
 * @desc This file is based on:
 * 	- the descriptions and Mathematica implementations of the field-polarised magnon model by M. Garst, 2016.
 * @license GPLv2 (see 'LICENSE' file)
 */

#include "fp.h"


// instantiation
template class FP<double, std::complex<double>>;
#ifdef __HACK_FULL_INST__
	template FP<double, std::complex<double>>::FP();
	template void FP<double, std::complex<double>>::SetG(double, double, double);
	template void FP<double, std::complex<double>>::SetCoords(double, double, double);
#endif



template<class t_real, class t_cplx>
FP<t_real, t_cplx>::FP()
{}


/**
 * set lattice vector and orthogonal projector
 */
template<class t_real, class t_cplx>
void FP<t_real, t_cplx>::SetG(t_real h, t_real k, t_real l)
{
	m_Grlu = tl2::make_vec<t_vec>({h,k,l});

	t_vec G = m_Grlu / tl2::veclen(m_Grlu);
	G = tl2::quat_vec_prod(m_rotCoord, G);
	m_projNeutron = tl2::unit_m<t_mat>(3) - tl2::outer<t_vec, t_mat>(G, G);
}


/**
 * rotates field to internal [001] convention
 */
template<class t_real, class t_cplx>
void FP<t_real, t_cplx>::SetCoords(t_real Bx, t_real By, t_real Bz)
{
	t_vec B = tl2::make_vec<t_vec>( {Bx, By, Bz} );
	t_quat quatB = tl2::rotation_quat(B, tl2::make_vec<t_vec>( {0, 0, 1} ));

	m_rotCoord = quatB;
}


/**
 * query the dispersion
 */
template<class t_real, class t_cplx>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
FP<t_real, t_cplx>::GetDisp(t_real h, t_real k, t_real l, t_real minE, t_real maxE) const
{
	constexpr auto imag = t_cplx(0,1);
	const auto ident2 = tl2::unit_m<t_mat_cplx>(2);
	const auto sigma = tl2::get_spin_matrices<ublas::matrix, std::vector, t_real>();
	const t_real dh0_shift = 0.007523;
	const t_real dh0 = dh0_shift + (m_B - m_Bc2) / m_Bc2;	// for given g_hoc


	t_vec Qvec = tl2::make_vec<t_vec>({ h, k, l });
	t_vec qvec = tl2::quat_vec_prod(m_rotCoord, Qvec) - tl2::quat_vec_prod(m_rotCoord, m_Grlu);
	qvec /= g_kh_rlu<t_real>(m_T);
	qvec = -qvec;

	t_real q=0., theta=0., phi=0.;
	std::tie(q, phi, theta) = tl2::cart_to_sph(qvec[0], qvec[1], qvec[2]);


	// eigensystems
	auto get_evecs = [&qvec, &q, &dh0, &imag, &sigma, &ident2](t_real phi) -> auto
	{
		// Hamiltonian
		t_cplx qp = qvec[0] + imag*qvec[1];
		t_cplx qm = qvec[0] - imag*qvec[1];

		t_mat_cplx H = (q*q + g_hoc<t_real>*q*q*q*q + 1. + dh0) * ident2;
		H += 2.*sigma[2] * qvec[2];
		H += g_chi<t_real>/(2.*q*q) * tl2::make_mat<t_mat_cplx>(
			{{ qm*qp, qm*qm },
			{ qp*qp, qm*qp }
		});


		// eigenvectors and -values
		std::vector<t_vec_cplx> evecs;
		std::vector<t_real> evals;
		bool ev_ok = tl2::eigenvecsel_herm<t_real>(tl2::prod_mm(sigma[2], H), evecs, evals, true);
		/*if(ev_ok)
		{
			t_cplx norm = tl2::mat_elem(evecs[0], sigma[2], evecs[1]);
			evecs[0] /= norm;
			evecs[1] /= norm;
		}*/

		return std::make_tuple(ev_ok, evals, evecs);
	};


	/*
	 * spectral weights
	 *
	 * transition probability to flip spin given by propagator, i.e. Green's function:
	 * G(x, x') = sum{j} |v_j'> <v_j| / E_j
	 * with eigenvectors |v_j> and eigenvalues E_j.
	 */
	auto get_weights = [&imag](const t_vec_cplx& evec_phi, const t_vec_cplx& evec_mphi, const t_mat_cplx& projNeutron, bool bProjNeutron) -> auto
	{
		t_mat_cplx kernel = tl2::outer<t_vec_cplx, t_mat_cplx>(evec_phi, evec_mphi);
		t_mat_cplx weight(3,3);

		const t_vec_cplx x = tl2::make_vec<t_vec_cplx>({ 1, 0, 0 });
		const t_vec_cplx y = tl2::make_vec<t_vec_cplx>({ 0, 1, 0 });

		for(int i=0; i<3; ++i)
		{
			t_vec_cplx veci = tl2::make_vec<t_vec_cplx>({ x[i]-imag*y[i], x[i]+imag*y[i] });

			for(int j=0; j<3; ++j)
			{
				t_vec_cplx vecj = tl2::make_vec<t_vec_cplx>({ x[j]-imag*y[j], x[j]+imag*y[j] });
				weight(i,j) = tl2::mat_elem(veci, kernel, vecj);
			}
		}


		// magnetic neutron scattering projections
		t_mat_cplx neutron_weight = weight;
		if(bProjNeutron)
		{
			neutron_weight = tl2::prod_mm(weight, projNeutron);
			neutron_weight = tl2::prod_mm(projNeutron, neutron_weight);
		}


		// polarisation projections
		t_mat_cplx matSF1 = tl2::prod_mm(get_chiralpol<t_mat_cplx>(1), neutron_weight);
		t_mat_cplx matSF2 = tl2::prod_mm(get_chiralpol<t_mat_cplx>(2), neutron_weight);
		t_mat_cplx matNSF = tl2::prod_mm(get_chiralpol<t_mat_cplx>(3), neutron_weight);

		constexpr const t_real weightScale = 0.5;
		t_cplx wAll = tl2::trace(neutron_weight) * weightScale;
		t_cplx wSF1 = tl2::trace(matSF1) * weightScale;
		t_cplx wSF2 = tl2::trace(matSF2) * weightScale;
		t_cplx wNSF = tl2::trace(matNSF) * weightScale;

		return std::make_tuple(wAll, wSF1, wSF2, wNSF);
	};


	auto [ok_phi, evals_phi, evecs_phi] = get_evecs(phi);
	auto [ok_mphi, evals_mphi, evecs_mphi] = get_evecs(-phi);
	if(!ok_phi || !ok_mphi)
	{
		std::vector<t_real> empty;
		return std::make_tuple(empty, empty, empty, empty, empty);
	}

	auto [wAll_p, wSF1_p, wSF2_p, wNSF_p] = get_weights(evecs_phi[0], evecs_mphi[0], m_projNeutron, m_bProjNeutron);
	auto [wAll_m, wSF1_m, wSF2_m, wNSF_m] = get_weights(-evecs_phi[1], evecs_mphi[1], m_projNeutron, m_bProjNeutron);


	std::vector<t_real> Es = { evals_phi[0]*g_g<t_real>*g_muB<t_real>*m_Bc2, evals_phi[1]*g_g<t_real>*g_muB<t_real>*m_Bc2 };
	std::vector<t_real> ws_all = { std::abs(wAll_p.real()), std::abs(wAll_m.real()) };
	std::vector<t_real> ws_SF1 = { std::abs(wSF1_p.real()), std::abs(wSF1_m.real()) };
	std::vector<t_real> ws_SF2 = { std::abs(wSF2_p.real()), std::abs(wSF2_m.real()) };
	std::vector<t_real> ws_NSF = { std::abs(wNSF_p.real()), std::abs(wNSF_m.real()) };

	return std::make_tuple(Es, ws_all, ws_SF1, ws_SF2, ws_NSF);
}
