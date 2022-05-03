/**
 * free energy and dynamics in the skyrmion phase
 * @author tweber@ill.fr
 * @date jul-18
 * @license GPLv2 (see 'LICENSE' file)
 */

// instantiation of skx module
#ifdef DEF_SKX_ORDER
	#pragma message("Skx Order: " __TL2_STRCONV(DEF_SKX_ORDER))
	template class Skx<double, std::complex<double>, DEF_SKX_ORDER>;

	#ifdef __HACK_FULL_INST__
		template Skx<double, std::complex<double>, DEF_SKX_ORDER>::Skx();
		template void Skx<double, std::complex<double>, DEF_SKX_ORDER>::SetG(double, double, double, bool);
		template void Skx<double, std::complex<double>, DEF_SKX_ORDER>::SetCoords(double, double, double, double, double, double);
	#endif
#endif


#ifndef SKX_USE_HOC
	#define SKX_USE_HOC 0
#endif
