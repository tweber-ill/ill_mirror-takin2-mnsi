/**
 * Helper functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date jul-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#ifndef __HELPER_H__
#define __HELPER_H__

#define USE_LAPACK
//#include "tlibs2/libs/math17.h"
//#include "tlibs2/libs/mag.h"
#include "tlibs2-extras/math17.h"
#include "tlibs2-extras/mag.h"
namespace ublas = boost::numeric::ublas;


/**
 * splits a vector into its x,y,z components
 */
template<class T = std::complex<double>>
std::tuple<T, T, T> split_vec3d(const ublas::vector<T>& vec)
{
	return std::make_tuple(vec[0], vec[1], vec[2]);
}



/**
 * array indexing including negative values
 */
template<class t_arr>
typename t_arr::value_type& get_comp(t_arr &arr, int idx)
{
	if(idx >= 0) return arr[idx];

	return arr[arr.size() + idx];
}



/**
 * matrix indexing including negative values
 */
template<class t_mat>
typename t_mat::value_type& get_comp(t_mat &mat, int idx1, int idx2)
{
	if(idx1 < 0) idx1 = int(mat.size1()) + idx1;
	if(idx2 < 0) idx2 = int(mat.size2()) + idx2;

	return mat(idx1, idx2);
}


/**
 * matrix indexing including negative values
 */
template<class t_mat>
const typename t_mat::value_type& get_comp(const t_mat &mat, int idx1, int idx2)
{
	if(idx1 < 0) idx1 = int(mat.size1()) + idx1;
	if(idx2 < 0) idx2 = int(mat.size2()) + idx2;

	return mat(idx1, idx2);
}


/**
 * 4-dim tensor indexing including negative values
 */
template<class t_arr>
typename t_arr::value_type& get_comp(t_arr& arr, int SIZE, int idx1, int idx2, int idx3, int idx4)
{
	if(idx1 < 0) idx1 = SIZE + idx1;
	if(idx2 < 0) idx2 = SIZE + idx2;
	if(idx3 < 0) idx3 = SIZE + idx3;
	if(idx4 < 0) idx4 = SIZE + idx4;

	return arr[idx1*SIZE*SIZE*SIZE +
		idx2*SIZE*SIZE +
		idx3*SIZE +
		idx4];
}


/**
 * 5-dim tensor indexing including negative values
 */
template<class t_arr>
typename t_arr::value_type& get_comp(t_arr& arr, int SIZE, int idx1, int idx2, int idx3, int idx4, int idx5)
{
	if(idx1 < 0) idx1 = SIZE + idx1;
	if(idx2 < 0) idx2 = SIZE + idx2;
	if(idx3 < 0) idx3 = SIZE + idx3;
	if(idx4 < 0) idx4 = SIZE + idx4;
	if(idx5 < 0) idx5 = SIZE + idx5;

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
typename t_arr::value_type& get_comp(t_arr& arr, int SIZE, int idx1, int idx2, int idx3, int idx4, int idx5, int idx6)
{
	if(idx1 < 0) idx1 = SIZE + idx1;
	if(idx2 < 0) idx2 = SIZE + idx2;
	if(idx3 < 0) idx3 = SIZE + idx3;
	if(idx4 < 0) idx4 = SIZE + idx4;
	if(idx5 < 0) idx5 = SIZE + idx5;
	if(idx6 < 0) idx6 = SIZE + idx6;

	return arr[idx1*SIZE*SIZE*SIZE*SIZE*SIZE +
		idx2*SIZE*SIZE*SIZE*SIZE +
		idx3*SIZE*SIZE*SIZE +
		idx4*SIZE*SIZE +
		idx5*SIZE +
		idx6];
}


template<class t_arr>
const typename t_arr::value_type& get_virt_comp(const t_arr& arr, int ORGSIZE, int VIRTSIZE, int ORDER, int idx1, int idx2, int idx3, int idx4)
{
	if(idx1 < 0) idx1 = VIRTSIZE + idx1;
	if(idx2 < 0) idx2 = VIRTSIZE + idx2;
	if(idx3 < 0) idx3 = VIRTSIZE + idx3;
	if(idx4 < 0) idx4 = VIRTSIZE + idx4;

	static const typename t_arr::value_type zero{};
	if((idx1>=ORDER && idx1<VIRTSIZE-ORDER-1) || (idx2>=ORDER && idx2<VIRTSIZE-ORDER-1) 
		|| (idx3>=ORDER && idx3<VIRTSIZE-ORDER-1) || (idx4>=ORDER && idx4<VIRTSIZE-ORDER-1))
		return zero;

	bool bLower1 = (idx1<ORDER);
	bool bLower2 = (idx2<ORDER);
	bool bLower3 = (idx3<ORDER);
	bool bLower4 = (idx4<ORDER);

	if(bLower1 && idx1 >= ORGSIZE) return zero;
	if(bLower2 && idx2 >= ORGSIZE) return zero;
	if(bLower3 && idx3 >= ORGSIZE) return zero;
	if(bLower4 && idx4 >= ORGSIZE) return zero;

	if(!bLower1)
	{
		idx1 -= VIRTSIZE-ORGSIZE;
		if(idx1 < 0) return zero;
	}
	if(!bLower2)
	{
		idx2 -= VIRTSIZE-ORGSIZE;
		if(idx2 < 0) return zero;
	}
	if(!bLower3)
	{
		idx3 -= VIRTSIZE-ORGSIZE;
		if(idx3 < 0) return zero;
	}
	if(!bLower4)
	{
		idx4 -= VIRTSIZE-ORGSIZE;
		if(idx4 < 0) return zero;
	}

	return arr[idx1*ORGSIZE*ORGSIZE*ORGSIZE + idx2*ORGSIZE*ORGSIZE + idx3*ORGSIZE + idx4];
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
	using t_cplx = typename t_mat::value_type;
	constexpr t_cplx j = t_cplx(0,1);

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
 * calculate energies and weights from Landau-Lifshitz Mcross and fluctuation matrices
 * Mcross rotates magnetisation
 *	- has 3 eigenvalues: 0, 1, -1:
 *		- 0: parallel magnetisation, long. fluctuation, not used
 *		- 1, -1: trans. fluctuations
 */
template<class t_mat_cplx, class t_vec_cplx, class t_cplx, class t_real>
std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
calc_weights(const t_mat_cplx& Mx, const t_mat_cplx& Fluc,
	bool bProjNeutron, const t_mat_cplx& projNeutron, const std::vector<t_mat_cplx>& polMat,
	t_real normfac=1, t_real E_scale_fac=1, t_real minE=-1, t_real maxE=-2,
	t_real eveps = 1e-6, t_real evlimit = 0.9995, t_real weighteps = 1e-6,
	bool bfilterzeroweight=0, int onlymode=-1, std::size_t MxsubMatRowBegin=0)
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
			Mx, Fluc, normfac, &eveps, &evlimit, 3, MxsubMatRowBegin, eveps);

		std::size_t numEVals = std::get<0>(eigs).size();
		std::vector<EW> EWs;
		EWs.reserve(numEVals);

		int iCurMode = 0;
		for(std::size_t ieval=0; ieval<numEVals; ++ieval)
		{
			auto E_weight = tl2::get_dynstrucfact_neutron<t_mat_cplx, t_vec_cplx, t_cplx, t_real>(
				std::get<0>(eigs)[ieval], std::get<1>(eigs)[ieval], std::get<2>(eigs)[ieval],
				bProjNeutron ? &projNeutron : nullptr, &polMat);

			// filter energies if requested
			t_real E_meV = std::get<0>(E_weight) * E_scale_fac;
			if(maxE >= minE && (E_meV < minE || E_meV > maxE))
					continue;

			// only consider one of the modes, if desired
			if(onlymode >=0 && onlymode != iCurMode++)
					continue;


			EW ew;
			ew.E = E_meV;
			ew.wUnpol = std::get<1>(E_weight)[0];
			ew.wSF1 = std::get<1>(E_weight)[1];
			ew.wSF2 = std::get<1>(E_weight)[2];
			ew.wNSF = std::get<1>(E_weight)[3];
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


#endif
