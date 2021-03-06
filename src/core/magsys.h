/**
 * Magnetic system interface
 * @author Tobias Weber <tweber@ill.fr>
 * @date jul-18
 * @license GPLv2 (see 'LICENSE' file)
 */

#ifndef __MAGSYS_H__
#define __MAGSYS_H__


#include <vector>
#include <complex>
#include <memory>
#include "helper.h"
#include "constants.h"


template<class t_real = double>
class HasTandB
{
public:
	virtual void SetB(t_real B) = 0;
	virtual void SetT(t_real T) = 0;
};


/**
 * interface for a magnetic system having a free energy and fourier components defined
 */
template<class t_real = double, class t_cplx = std::complex<t_real>, int ORDER_FOURIER=2>
class MagSystem : public HasTandB<t_real>
{
public:
	using value_type = t_real;
	using value_type_c = t_cplx;

public:
	virtual t_real F() = 0;

	virtual void SetFourier(const std::vector<ublas::vector<t_cplx>> &fourier) = 0;
	virtual const std::vector<ublas::vector<t_cplx>> &GetFourier() const = 0;

	virtual bool minimise(int iMaxOrder = 9999,
		bool bFixXR=0, bool bFixYR=0, bool bFixZR=0, bool bFixXI=0, bool bFixYI=0, bool bFixZI=0);

	virtual bool SaveStates(const char *file, int iMaxOrder = 9999,
		bool bFixXR=0, bool bFixYR=0, bool bFixZR=0, bool bFixXI=0, bool bFixYI=0, bool bFixZI=0) const;

	void SetDebug(bool b) { m_debug = b; }

	virtual std::shared_ptr<MagSystem<t_real, t_cplx, ORDER_FOURIER>> copyCastSys() const = 0;

protected:
	bool m_debug = false;
};


/**
 * interface for a magnon dynamics
 */
template<class t_real = double, class t_cplx = std::complex<t_real>>
class MagDynamics : public HasTandB<t_real>
{
public:
	// unpol, SF-+, SF+-, NSF
	virtual std::tuple<std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>, std::vector<t_real>>
		GetDisp(t_real h, t_real k, t_real l, t_real minE=-1., t_real maxE=-2.) const = 0;

	virtual std::shared_ptr<MagDynamics<t_real, t_cplx>> copyCastDyn() const = 0;
};


#endif
