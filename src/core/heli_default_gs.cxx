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
	{ 0, 0 }, { 0, -0 }, { 11.154587, 0 },
	{ 6.8518185, 6.8518554 }, { 6.8518185, -6.8518554 }, { 0.0074620428, 0 },
	{ -0.0027917535, -0.0031683066 }, { -0.0027917535, 0.0031683066 }, { 0.00094754919, 0 },
	{ -0.00056589914, -0.00011888157 }, { -0.00056589914, 0.00011888157 }, { 9.8035831e-05, 0 },
	{ 0.00012901039, -0.0005550897 }, { 0.00012901039, 0.0005550897 }, { 0.0014994666, 0 },
	{ -0.00029251049, -0.0004032893 }, { -0.00029251049, 0.0004032893 }, { -0.0015262806, 0 },
	{ 0.00061386564, -0.00015588607 }, { 0.00061386564, 0.00015588607 }, { -0.00082875038, 0 },
	{ 0.00035018072, 0.00017051901 }, { 0.00035018072, -0.00017051901 }, { -0.0027675351, 0 },
	{ 0.00030682248, 0.00098032048 }, { 0.00030682248, -0.00098032048 }, { -0.0014957615, 0 },
	{ 0.00046924383, 0.00046777349 }, { 0.00046924383, -0.00046777349 }, { 0.00090007733, 0 },
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

