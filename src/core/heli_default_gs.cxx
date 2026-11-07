/**
 * default helimagnetic ground state
 * @author tweber@ill.fr
 * @date apr-22
 * @license GPLv2 (see 'LICENSE' file)
 */


/**
 * temperature and field (both in theoretical units),
 * for which the ground state has been calculated
 */
static const t_real _heligs_T = -1000.;
static const t_real _heligs_B = 23.9283;


/**
 * output from heli_gs.cpp
 */
static const std::vector<std::complex<t_real>> _heligs_allcomps =
{{
	{ 0, 0 }, { 0, -0 }, { 11.874347, 0 }, // 0
	{ 6.6989227, 6.6989043 }, { 6.6989227, -6.6989043 }, { -1.3652825e-05, 0 }, // 1
	{ -7.6951827e-06, 1.9545602e-05 }, { -7.6951827e-06, -1.9545602e-05 }, { -4.9845865e-06, 0 }, // 2
	{ 1.0236135e-05, -2.4457819e-05 }, { 1.0236135e-05, 2.4457819e-05 }, { 1.7881011e-05, 0 }, // 3
	{ -9.290194e-06, 2.0686519e-05 }, { -9.290194e-06, -2.0686519e-05 }, { -2.8708125e-05, 0 }, // 4
	{ 3.189548e-06, -1.0259606e-05 }, { 3.189548e-06, 1.0259606e-05 }, { 1.5167963e-05, 0 }, // 5
	{ -2.5547014e-06, 1.8488461e-05 }, { -2.5547014e-06, -1.8488461e-05 }, { -8.4996825e-06, 0 }, // 6
	{ -9.2606592e-06, 4.9391666e-06 }, { -9.2606592e-06, -4.9391666e-06 }, { -4.8434831e-06, 0 }, // 7
	{ 1.3250869e-05, 1.1433526e-05 }, { 1.3250869e-05, -1.1433526e-05 }, { 1.3447745e-05, 0 }, // 8
	{ 1.3250869e-05, -1.1433526e-05 }, { 1.3250869e-05, 1.1433526e-05 }, { 1.3447745e-05, -0 }, // -8
	{ -9.2606592e-06, -4.9391666e-06 }, { -9.2606592e-06, 4.9391666e-06 }, { -4.8434831e-06, -0 }, // -7
	{ -2.5547014e-06, -1.8488461e-05 }, { -2.5547014e-06, 1.8488461e-05 }, { -8.4996825e-06, -0 }, // -6
	{ 3.189548e-06, 1.0259606e-05 }, { 3.189548e-06, -1.0259606e-05 }, { 1.5167963e-05, -0 }, // -5
	{ -9.290194e-06, -2.0686519e-05 }, { -9.290194e-06, 2.0686519e-05 }, { -2.8708125e-05, -0 }, // -4
	{ 1.0236135e-05, 2.4457819e-05 }, { 1.0236135e-05, -2.4457819e-05 }, { 1.7881011e-05, -0 }, // -3
	{ -7.6951827e-06, -1.9545602e-05 }, { -7.6951827e-06, 1.9545602e-05 }, { -4.9845865e-06, -0 }, // -2
	{ 6.6989227, -6.6989043 }, { 6.6989227, 6.6989043 }, { -1.3652825e-05, -0 }, // -1
}};


template<class t_vec_cplx>
std::tuple<
	t_real,                  // T
	t_real,                  // B
	std::vector<t_vec_cplx>  // fourier components
>
_get_heli_gs()
{
	// set initial fourier components
	std::vector<t_vec_cplx> fourier;
	fourier.reserve(_heligs_allcomps.size()/3);

	for(std::size_t comp=0; comp<_heligs_allcomps.size(); comp+=3)
	{
		fourier.emplace_back(tl2::make_vec<t_vec_cplx>(
		{
			_heligs_allcomps[comp],
			_heligs_allcomps[comp+1],
			_heligs_allcomps[comp+2]
		}));
	}

	return std::make_tuple(_heligs_T, _heligs_B, fourier);
}
