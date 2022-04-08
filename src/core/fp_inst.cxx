/**
 * dynamics in the field-polarised phase
 * @author Tobias Weber <tweber@ill.fr>
 * @date dec-16, sep-18
 * @license GPLv2 (see 'LICENSE' file)
 */

// instantiation of fp module
template class FP<double, std::complex<double>>;
#ifdef __HACK_FULL_INST__
	template FP<double, std::complex<double>>::FP();
	template void FP<double, std::complex<double>>::SetG(double, double, double);
	template void FP<double, std::complex<double>>::SetCoords(double, double, double, double, double, double);
#endif
