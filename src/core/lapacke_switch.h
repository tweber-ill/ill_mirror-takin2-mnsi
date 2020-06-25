/**
 * Lapacke chooser
 * @author tweber@ill.fr
 * @date sep-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#ifndef __LAPACKE_SWITCH_H__
#define __LAPACKE_SWITCH_H__

#if __has_include(<mkl_lapacke.h>)
	#include <mkl_lapacke.h>
	#pragma message "Using MKL Lapacke."
#else
	#include <lapacke.h>
	#pragma message "Using Vanilla Lapacke."
#endif

#endif
