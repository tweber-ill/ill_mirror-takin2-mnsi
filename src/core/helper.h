/**
 * helper functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date jul-18
 * @license GPLv2 (see 'LICENSE' file)
 * @desc This file implements the theoretical skyrmion model by M. Garst and J. Waizner, references:
 *	- M. Garst, J. Waizner, and D. Grundler, J. Phys. D: Appl. Phys. 50, 293002 (2017), https://doi.org/10.1088/1361-6463/aa7573
 *	- J. Waizner, PhD thesis (2017), Universität zu Köln, https://kups.ub.uni-koeln.de/7937/
 *	- Personal communications with M. Garst, 2017-2020.
 */

#ifndef __HELPER_H__
#define __HELPER_H__

#define USE_LAPACK
#include "tlibs2/libs/math17.h"
#include "tlibs2/libs/mag.h"
namespace ublas = boost::numeric::ublas;

#include <type_traits>


/**
 * splits a vector into its x,y,z components
 */
template<class T = std::complex<double>>
std::tuple<T, T, T> split_vec3d(const ublas::vector<T>& vec)
{
	return std::make_tuple(vec[0], vec[1], vec[2]);
}


/**
 * index sequence including negative indices: 0, 1, ..., ORDER, -ORDER, -ORDER+1, ..., -1
 */
template<class t_int = int>
t_int abs_to_rel_idx(t_int h_idx, t_int ORDER)
{
	const t_int SIZE = 2*ORDER + 1;
	return (h_idx <= ORDER) ? h_idx : h_idx-SIZE;
}


/**
 * array indexing including negative indices
 */
template<class t_int = int>
t_int rel_to_abs_idx(t_int h, t_int SIZE)
{
	return (h >= 0) ? h : SIZE + h;
}


/**
 * array indexing including negative values
 */
template<class t_arr>
std::conditional_t<std::is_const_v<t_arr>,
	const typename t_arr::value_type&,
	typename t_arr::value_type&>
get_comp(t_arr &arr, int idx)
{
	idx = rel_to_abs_idx(idx, static_cast<int>(arr.size()));
	return arr[idx];
}


/**
 * matrix indexing including negative values
 */
template<class t_mat>
std::conditional_t<std::is_const_v<t_mat>,
	const typename t_mat::value_type&,
	typename t_mat::value_type&>
get_comp(t_mat &mat, int idx1, int idx2)
{
	idx1 = rel_to_abs_idx(idx1, static_cast<int>(mat.size1()));
	idx2 = rel_to_abs_idx(idx2, static_cast<int>(mat.size2()));

	return mat(idx1, idx2);
}


/**
 * 2-dim tensor indexing including negative values
 */
template<class t_arr>
typename t_arr::value_type& get_comp(t_arr& arr, int SIZE,
	int idx1, int idx2)
{
	idx1 = rel_to_abs_idx(idx1, SIZE);
	idx2 = rel_to_abs_idx(idx2, SIZE);

	return arr[idx1*SIZE + idx2];
}


/**
 * 4-dim tensor indexing including negative values
 */
template<class t_arr>
typename t_arr::value_type& get_comp(t_arr& arr, int SIZE,
	int idx1, int idx2, int idx3, int idx4)
{
	idx1 = rel_to_abs_idx(idx1, SIZE);
	idx2 = rel_to_abs_idx(idx2, SIZE);
	idx3 = rel_to_abs_idx(idx3, SIZE);
	idx4 = rel_to_abs_idx(idx4, SIZE);

	return arr[idx1*SIZE*SIZE*SIZE +
		idx2*SIZE*SIZE +
		idx3*SIZE +
		idx4];
}


/**
 * 5-dim tensor indexing including negative values
 */
template<class t_arr>
typename t_arr::value_type& get_comp(t_arr& arr, int SIZE,
	int idx1, int idx2, int idx3, int idx4, int idx5)
{
	idx1 = rel_to_abs_idx(idx1, SIZE);
	idx2 = rel_to_abs_idx(idx2, SIZE);
	idx3 = rel_to_abs_idx(idx3, SIZE);
	idx4 = rel_to_abs_idx(idx4, SIZE);
	idx5 = rel_to_abs_idx(idx5, SIZE);

	return arr[idx1*SIZE*SIZE*SIZE*SIZE +
		idx2*SIZE*SIZE*SIZE +
		idx3*SIZE*SIZE +
		idx4*SIZE +
		idx5];
}


/**
 * 6-dim tensor indexing including negative values
 */
template<class t_arr>
typename t_arr::value_type& get_comp(t_arr& arr, int SIZE,
	int idx1, int idx2, int idx3, int idx4, int idx5, int idx6)
{
	idx1 = rel_to_abs_idx(idx1, SIZE);
	idx2 = rel_to_abs_idx(idx2, SIZE);
	idx3 = rel_to_abs_idx(idx3, SIZE);
	idx4 = rel_to_abs_idx(idx4, SIZE);
	idx5 = rel_to_abs_idx(idx5, SIZE);
	idx6 = rel_to_abs_idx(idx6, SIZE);

	return arr[idx1*SIZE*SIZE*SIZE*SIZE*SIZE +
		idx2*SIZE*SIZE*SIZE*SIZE +
		idx3*SIZE*SIZE*SIZE +
		idx4*SIZE*SIZE +
		idx5*SIZE +
		idx6];
}


/**
 * index into extended system
 */
template<class t_arr>
const typename t_arr::value_type& get_ext_comp(bool use_ext_sys,
	t_arr& arr, int ARRSIZE, int ORDER, int MAXORDER,
	int idx1, int idx2, int idx3, int idx4)
{
	if(!use_ext_sys)
		return get_comp(arr, ARRSIZE, idx1, idx2, idx3, idx4);

	const int MAXSIZE = 2*MAXORDER + 1;

	// negative indices
	idx1 = rel_to_abs_idx(idx1, MAXSIZE);
	idx2 = rel_to_abs_idx(idx2, MAXSIZE);
	idx3 = rel_to_abs_idx(idx3, MAXSIZE);
	idx4 = rel_to_abs_idx(idx4, MAXSIZE);

	static const typename t_arr::value_type zero{};
	const int diffsize = MAXSIZE-ARRSIZE;

	if(idx1 >= ORDER && (idx1 < MAXSIZE-ORDER-1 || idx1 < diffsize)) return zero;
	if(idx2 >= ORDER && (idx2 < MAXSIZE-ORDER-1 || idx2 < diffsize)) return zero;
	if(idx3 >= ORDER && (idx3 < MAXSIZE-ORDER-1 || idx3 < diffsize)) return zero;
	if(idx4 >= ORDER && (idx4 < MAXSIZE-ORDER-1 || idx4 < diffsize)) return zero;

	if(idx1 < ORDER && idx1 >= ARRSIZE) return zero;
	if(idx2 < ORDER && idx2 >= ARRSIZE) return zero;
	if(idx3 < ORDER && idx3 >= ARRSIZE) return zero;
	if(idx4 < ORDER && idx4 >= ARRSIZE) return zero;

	if(idx1 >= ORDER) idx1 -= diffsize;
	if(idx2 >= ORDER) idx2 -= diffsize;
	if(idx3 >= ORDER) idx3 -= diffsize;
	if(idx4 >= ORDER) idx4 -= diffsize;

	return arr[idx1*ARRSIZE*ARRSIZE*ARRSIZE +
		idx2*ARRSIZE*ARRSIZE +
		idx3*ARRSIZE +
		idx4];
}


/**
 * get chiral basis vectors
 */
template<class t_vec = ublas::vector<std::complex<double>>>
t_vec get_chiralbasis(int which)
{
	using t_cplx = typename t_vec::value_type;
	constexpr t_cplx j = t_cplx(0,1);

	t_cplx n = std::sqrt<typename t_cplx::value_type>(2.);

	if(which == 1)			// SF-+
		return tl2::make_vec<t_vec>( {1, +j, 0} ) / n;
	else if(which == 2)		// SF+-
		return tl2::make_vec<t_vec>( {1, -j, 0} ) / n;
	else if(which == 3)		// NSF
		return tl2::make_vec<t_vec>( {0, 0, 1} );

	// none
	return tl2::make_vec<t_vec>( {0, 0, 0} );
}


/**
 * get chiral basis matrix
 */
template<class t_mat = ublas::matrix<std::complex<double>>,
	class t_vec = ublas::vector<std::complex<double>>>
t_mat get_chiralbasismat()
{
	t_vec vec1 = get_chiralbasis<t_vec>(1);
	t_vec vec2 = get_chiralbasis<t_vec>(2);
	t_vec vec3 = get_chiralbasis<t_vec>(3);

	return tl2::row_matrix<t_mat, t_vec>({vec1, vec2, vec3});
}


/**
 * get polarisation matrix
 */
template<class t_mat = ublas::matrix<std::complex<double>>>
t_mat get_polmat(int which)
{
	// all channels
	if(which < 1 || which > 3) return tl2::unit_m<t_mat>(3);

	t_mat mat = tl2::zero_m<t_mat>(3, 3);
	mat(which-1, which-1) = 1;	// SF-+, SF+-, NSF

	return mat;
}


/**
 * get polarisation matrix in chiral basis
 */
template<class t_mat = ublas::matrix<std::complex<double>>>
t_mat get_chiralpol(int which)
{
	using t_cplx = typename t_mat::value_type;
	constexpr t_cplx j = t_cplx(0,1);

	if(which == 1)		// SF-+
		return tl2::make_mat<t_mat>({{0.5, +0.5*j, 0}, {-0.5*j, 0.5, 0}, {0, 0, 0}});
	else if(which == 2)	// SF+-
		return tl2::make_mat<t_mat>({{0.5, -0.5*j, 0}, {+0.5*j, 0.5, 0}, {0, 0, 0}});
	else if(which == 3)	// NSF
		return tl2::make_mat<t_mat>({{0, 0, 0}, {0, 0, 0}, {0, 0, 1}});

	// all channels
	return tl2::unit_m<t_mat>(3);
}


/**
 * calculate energies and weights from Landau-Lifshitz M-cross and fluctuation matrices
 */
template<class t_mat_cplx, class t_vec_cplx, class t_cplx, class t_real>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
calc_weights(const t_mat_cplx& Mx, const t_mat_cplx& Fluc,
	const t_mat_cplx& projNeutron, const std::vector<t_mat_cplx>& polMat,
	t_real E_scale_fac = 1., t_real w_scale_fac = 1., t_real minE = -1, t_real maxE = -2,
	t_real eveps = 1e-6, t_real evlimit = 0.9995, t_real weighteps = 1e-6,
	bool bfilterzeroweight = 0, int onlymode = -1, std::size_t MxsubMatRowBegin = 0)
{
	try
	{
		// energies and spectral weights
		struct EW
		{
			t_real E;
			t_real wUnpol;
			t_real wSF1, wSF2, wNSF;
		};

		auto eigs = tl2::calc_dynstrucfact_landau<t_mat_cplx, t_vec_cplx, t_cplx, t_real>(
			Mx, Fluc, &eveps, &evlimit, 3, MxsubMatRowBegin, eveps);

		std::size_t numEVals = std::get<0>(eigs).size();
		std::vector<EW> EWs;
		EWs.reserve(numEVals);

		bool bProjNeutron = true;
		int iCurMode = 0;
		for(std::size_t ieval=0; ieval<numEVals; ++ieval)
		{
			auto E_weight = tl2::get_dynstrucfact_neutron<t_mat_cplx, t_vec_cplx, t_cplx, t_real>(
				std::get<0>(eigs)[ieval],  std::get<2>(eigs)[ieval],
				bProjNeutron ? &projNeutron : nullptr, &polMat);

			t_real E_meV = std::get<0>(E_weight) * E_scale_fac;

			// filter energies if requested
			if(maxE >= minE && (E_meV < minE || E_meV > maxE))
				continue;

			// only consider one of the modes, if desired
			if(onlymode >=0 && onlymode != iCurMode++)
				continue;

			EW ew;
			ew.E = E_meV;
			ew.wUnpol = std::get<1>(E_weight)[0] * w_scale_fac;
			ew.wSF1 = std::get<1>(E_weight)[1] * w_scale_fac;
			ew.wSF2 = std::get<1>(E_weight)[2] * w_scale_fac;
			ew.wNSF = std::get<1>(E_weight)[3] * w_scale_fac;
			EWs.emplace_back(ew);
		}


		std::stable_sort(EWs.begin(), EWs.end(), [](const EW& ew1, const EW& ew2) -> bool
		{ return ew1.E < ew2.E; });

		std::vector<t_real> Es, wsUnpol, wsSF1, wsSF2, wsNSF;
		for(const EW& ew : EWs)
		{
			if(bfilterzeroweight && ew.wUnpol < weighteps)
				continue;

			Es.push_back(ew.E);
			wsUnpol.push_back(ew.wUnpol);
			wsSF1.push_back(ew.wSF1);
			wsSF2.push_back(ew.wSF2);
			wsNSF.push_back(ew.wNSF);
		}

		return std::make_tuple(Es, wsUnpol, wsSF1, wsSF2, wsNSF);
	}
	catch(const std::exception& ex)
	{
		tl2::log_err(ex.what());

		static const std::vector<t_real> empty;
		return std::make_tuple(empty, empty, empty, empty, empty);
	}
}


/**
 * appends multiple containers to one another
 */
template<class t_cont, std::size_t... seq>
void insert_vals(t_cont& val0, const t_cont& val1, std::index_sequence<seq...>)
{
	( std::get<seq>(val0).insert(std::get<seq>(val0).end(),
		std::get<seq>(val1).begin(), std::get<seq>(val1).end()), ... );
}


template<class t_mat1, class t_mat2>
void assign_or_add(t_mat1& matDst, const t_mat2& mat)
{
	if(!matDst.size1())
		matDst = mat;
	else
		matDst += mat;
}


/**
 * move the momentum transfer away from the lattice vector to avoid divergences
 */
template<class t_real>
void avoid_G(t_real& qh, t_real& qk, t_real& ql, t_real eps)
{
	if(tl2::float_equal<t_real>(qh, 0., eps) &&
		tl2::float_equal<t_real>(qk, 0., eps) &&
		tl2::float_equal<t_real>(ql, 0., eps))
	{
		qh += eps;
		qk += eps;
		ql += eps;
	}
}


/**
 * get the next lattice miller index
 */
template<class t_int = int, class t_real = double>
t_int lattidx(t_real h)
{
	return static_cast<t_int>(std::round(h));
}


#endif
