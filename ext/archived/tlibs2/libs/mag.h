/**
 * tlibs2 -- magnon helper functions moved from the main tlibs2 repo: https://code.ill.fr/scientific-software/takin/tlibs2
 * @author Tobias Weber <tweber@ill.fr>
 * @date 30-may-2020
 * @license GPLv2 or GPLv3, see 'LICENSE' file
 * @desc tlibs forked on 7-Nov-2018 from the privately developed "tlibs" project (https://github.com/t-weber/tlibs).
 *
 * @desc This file is based on the theoretical helimagnon and skyrmion models by M. Garst and J. Waizner, see:
 *	- personal communications with M. Garst, 2017-2020.
 *	- M. Garst, J. Waizner, and D. Grundler, J. Phys. D: Appl. Phys. 50, 293002 (2017), https://doi.org/10.1088/1361-6463/aa7573
 *	- J. Waizner, PhD thesis (2017), Universität zu Köln, https://kups.ub.uni-koeln.de/7937/
 *	- M. Kugler, G. Brandl, et al., Phys. Rev. Lett. 115, 097203 (2015),https://doi.org/10.1103/PhysRevLett.115.097203
 *	- T Schwarze, J. Waizner, et al., Nat. Mater. 14, pp. 478–483 (2015), https://doi.org/10.1038/nmat4223
 *
 * ----------------------------------------------------------------------------
 * This file is a (former) part of tlibs (DOI: 10.5281/zenodo.5717779).
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 or version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_PHYS_MAG__
#define __TLIBS2_PHYS_MAG__

#include "math17.h"


// default optimisation flags
#ifndef MX_IS_HERM
	#define MX_IS_HERM 1
#endif
#ifndef MXX_IS_DIAG
	#define MXX_IS_DIAG 1
#endif
#ifndef INTERACTMAT_IS_HERM
	// warning: setting this to "1" leads to false low-q eigenvectors!
	#define INTERACTMAT_IS_HERM 0
#endif


#if !defined(__TL2_STRCONV0) && !defined(__TL2_STRCONV)
	#define __TL2_STRCONV0(__DEF) #__DEF
	#define __TL2_STRCONV(__DEF) __TL2_STRCONV0(__DEF)
#endif

#pragma message("MX_IS_HERM: " __TL2_STRCONV(MX_IS_HERM))
#pragma message("MXX_IS_DIAG: " __TL2_STRCONV(MXX_IS_DIAG))
#pragma message("INTERACTMAT_IS_HERM: " __TL2_STRCONV(INTERACTMAT_IS_HERM))


namespace tl2 {

/**
 * Calculates energies and dynamical structure factors from Landau-Lifshitz (M x) and fluctuation matrices.
 * Implements the mathematical formalism by M. Garst, J. Waizner, et al., references:
 *	- personal communications with M. Garst, 2017-2020.
 *	- https://doi.org/10.1088/1361-6463/aa7573
 *	- https://doi.org/10.1038/nmat4223 (supplement)
 *	- https://kups.ub.uni-koeln.de/7937/
 */
template<class t_mat_cplx, class t_vec_cplx, class t_cplx, class t_real>
std::tuple<std::vector<t_cplx>, std::vector<t_vec_cplx>, std::vector<t_mat_cplx>>
calc_dynstrucfact_landau(const t_mat_cplx& Mx, const t_mat_cplx& Fluc,
	const t_real* mineval = nullptr, const t_real *maxeval = nullptr,
	std::size_t MxsubMatSize = 3, std::size_t MxsubMatRowBegin = 0, t_real eps = 1e-6)
{
	constexpr t_cplx imag = t_cplx(0,1);

	// calculate Mx eigenvalues
	std::vector<t_vec_cplx> Mxevecs;
	std::vector<t_cplx> Mxevals;
	{
#if MX_IS_HERM == 1
		std::vector<t_real> _Mxevals;
		if(!eigenvecsel_herm<t_real>(-imag*Mx, Mxevecs, _Mxevals, true, -1., -2., eps))
			throw std::runtime_error("Mx eigenvector determination failed!");
		for(t_real d : _Mxevals)
			Mxevals.emplace_back(t_cplx(0., d));
#else
		if(!eigenvec_cplx<t_real>(Mx, Mxevecs, Mxevals, true))
			throw std::runtime_error("Mx eigenvector determination failed!");
#endif

		// filter eigenvalues
		auto maxelem = std::max_element(Mxevals.begin(), Mxevals.end(),
			[](const t_cplx& x, const t_cplx& y) -> bool
			{ return std::abs(x.imag()) < std::abs(y.imag()); });

		std::vector<t_vec_cplx> Mxevecs_new;
		std::vector<t_cplx> Mxevals_new;

		for(std::size_t elem=0; elem<Mxevals.size(); ++elem)
		{
			// upper eigenvalue limit
			if(maxeval && std::abs(Mxevals[elem].imag()) < std::abs(*maxelem)**maxeval)
				continue;
			// lower eigenvalue limit
			if(mineval && std::abs(Mxevals[elem].imag()) < *mineval)
				continue;
			Mxevecs_new.push_back(Mxevecs[elem]);
			Mxevals_new.push_back(Mxevals[elem]);
		}

		Mxevecs = std::move(Mxevecs_new);
		Mxevals = std::move(Mxevals_new);
	}


	// convert to eigenvector matrix
	t_mat_cplx MxEvecs(Mx.size1(), Mxevecs.size());
	for(std::size_t idx1=0; idx1<MxEvecs.size1(); ++idx1)
		for(std::size_t idx2=0; idx2<MxEvecs.size2(); ++idx2)
			MxEvecs(idx1, idx2) = Mxevecs[idx2][idx1];

	t_mat_cplx MxEvecsH = hermitian(MxEvecs);
	t_mat_cplx MxEvecs3 = submatrix_wnd<t_mat_cplx>(MxEvecs, MxsubMatSize, Mxevecs.size(), MxsubMatRowBegin, 0);
	t_mat_cplx MxEvecsH3 = hermitian(MxEvecs3);

	// transform fluctuation matrix into Mx eigenvector system
	t_mat_cplx invsuscept = prod_mm(Fluc, MxEvecs);
	invsuscept = prod_mm(MxEvecsH, invsuscept);

	// transform Mx into Mx eigenvector system
	// Mxx is diagonal with this construction => directly use Mxevals
#if MXX_IS_DIAG == 1
	t_mat_cplx Mxx = imag * diag_matrix<t_mat_cplx>(Mxevals);
#else
	t_mat_cplx Mxx = prod_mm(Mx, MxEvecs);
	Mxx = imag * prod_mm(MxEvecsH, Mxx);
#endif

	// Landau-Lifshitz: d/dt dM = -Mx B_mean, B_mean = -chi^(-1) * dM
	// E = EVals{ i Mx chi^(-1) }
	// chi_dyn^(-1) = i*E*Mx^(-1) + chi^(-1)
	// Mx*chi_dyn^(-1) = i*E + Mx*chi^(-1)
	t_mat_cplx Interactmat = prod_mm(Mxx, invsuscept);
	std::vector<t_vec_cplx> Interactevecs;
	std::vector<t_cplx> Interactevals;
#if INTERACTMAT_IS_HERM == 1
	std::vector<t_real> _Interactevals;
	if(!eigenvecsel_herm<t_real>(Interactmat, Interactevecs, _Interactevals, true, -1., -2., eps))
		throw std::runtime_error("Mxx eigenvector determination failed!");
	for(t_real d : _Interactevals)
		Interactevals.emplace_back(t_cplx(d, 0.));
#else
	if(!eigenvec_cplx<t_real>(Interactmat, Interactevecs, Interactevals, true))
		throw std::runtime_error("Mxx eigenvector determination failed!");
#endif

	std::vector<t_mat_cplx> Interactemats;
	Interactemats.reserve(Interactevals.size());

	for(std::size_t iInteract=0; iInteract<Interactevals.size(); ++iInteract)
	{
		const t_vec_cplx& evec = Interactevecs[iInteract];
		const t_vec_cplx evec_scale = prod_mv(Mxx, evec);

		auto matOuter = outer_cplx<t_vec_cplx, t_mat_cplx>(evec, evec);
		matOuter /= inner_cplx<t_vec_cplx>(evec, evec_scale);

		t_mat_cplx emat = prod_mm(matOuter, MxEvecsH3);
		emat = prod_mm(MxEvecs3, emat);
		Interactemats.emplace_back(std::move(emat));
	}

	return std::make_tuple(Interactevals, Interactevecs, Interactemats);
}



/**
 * Gets the dynamical structure factors from the eigenvectors calculated using calc_dynstrucfact_landau.
 * Implements the mathematical formalism by M. Garst, J. Waizner, et al., references:
 *	- personal communications with M. Garst, 2017-2020.
 *	- https://doi.org/10.1088/1361-6463/aa7573
 *	- https://doi.org/10.1038/nmat4223 (supplement)
 *	- https://kups.ub.uni-koeln.de/7937/
 */
template<class t_mat_cplx, class t_vec_cplx, class t_cplx, class t_real>
std::tuple<t_real, std::vector<t_real>>
get_dynstrucfact_neutron(
	const t_cplx& eval, const t_mat_cplx& _emat,
	const t_mat_cplx* projNeutron = nullptr,
	const std::vector<t_mat_cplx>* pol = nullptr)
{
	t_real E = eval.real();
	t_mat_cplx emat = _emat;

	// magnetic neutron scattering orthogonal projector: projNeutron = 1 - |G><G|
	if(projNeutron)
	{
		emat = prod_mm(emat, *projNeutron);
		emat = prod_mm(*projNeutron, emat);
	}

	std::vector<t_real> sfacts;

	// unpolarised structure factor
	sfacts.push_back(std::abs(trace(emat).real()));

	// polarised structure factors
	if(pol)
	{
		for(const t_mat_cplx& polmat : *pol)
		{
			t_mat_cplx matSF = prod_mm(polmat, emat);
			sfacts.push_back(std::abs(trace(matSF).real()));
		}
	}

	return std::make_tuple(E, sfacts);
}

}
#endif
