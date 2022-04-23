/**
 * default helimagnetic ground state
 * @author tweber@ill.fr
 * @date apr-22
 * @license GPLv2 (see 'LICENSE' file)
 */


/**
 * output from heli_gs.cpp
 */
static const std::vector<std::complex<t_real>> _heligs_allcomps =
{{
	{ 0, 0 }, { 0, -0 }, { 11.185589, 0 },
	{ 6.844417, 6.8466441 }, { 6.844417, -6.8466441 }, { 0.0059896075, 0 },
	{ -0.0023120373, -0.0026330922 }, { -0.0023120373, 0.0026330922 }, { 0.0018642475, 0 },
	{ 0.00013396703, 0.00067249845 }, { 0.00013396703, -0.00067249845 }, { -0.0005280514, 0 },
	{ 0.00031996453, -0.00026956086 }, { 0.00031996453, 0.00026956086 }, { 0.0024771008, 0 },
	{ -0.00088895265, -0.00056444692 }, { -0.00088895265, 0.00056444692 }, { -0.0014874747, 0 },
	{ 0.00067106976, -0.00015224015 }, { 0.00067106976, 0.00015224015 }, { 0.0014053083, 0 },
	{ -0.00039470042, -0.0003362198 }, { -0.00039470042, 0.0003362198 }, { -0.0034247079, 0 },
	{ 0.00058321402, 0.0011575263 }, { 0.00058321402, -0.0011575263 }, { 0.00016320065, 0 },
	{ -8.9200078e-06, -2.6656354e-05 }, { -8.9200078e-06, 2.6656354e-05 }, { 0.00080280193, 0 },
}};


template<class t_vec_cplx>
std::vector<t_vec_cplx> _get_heli_gs()
{
	using t_cplx = typename t_vec_cplx::value_type;
	constexpr auto imag = t_cplx(0, 1);

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

	return fourier;
}

