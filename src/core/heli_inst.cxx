/**
 * free energy and dynamics in the helimagnetic phase
 * @author tweber@ill.fr
 * @date mid-16, jul-18
 * @license GPLv2 (see 'LICENSE' file)
 */

// instantiation of heli module
#ifdef DEF_HELI_ORDER
	#pragma message("Heli Order: " __TL2_STRCONV(DEF_HELI_ORDER))
	template class Heli<double, std::complex<double>, DEF_HELI_ORDER>;

	#ifdef __HACK_FULL_INST__
		template Heli<double, std::complex<double>, DEF_HELI_ORDER>::Heli();
		template void Heli<double, std::complex<double>, DEF_HELI_ORDER>::SetCoords(double, double, double, double, double, double);
	#endif
#endif

#ifdef HELI_USE_HOC
	#pragma message("Heli HOC: " __TL2_STRCONV(HELI_USE_HOC))
#endif
