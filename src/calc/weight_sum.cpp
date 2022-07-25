/**
 * calculates the integrated weights/energies for setup 3
 * @author Tobias Weber <tweber@ill.fr>
 * @date jun-20
 * @license GPLv2 (see 'LICENSE' file)
 *
 * TODO: take into account the entire sector, not just the central path
 */

#include "core/heli.h"
#include "core/skx.h"
#include "tlibs2/libs/phys.h"

#include <fstream>

#include <boost/histogram.hpp>
namespace hist = boost::histogram;

using t_real = double;
using t_cplx = std::complex<t_real>;
using t_vec = ublas::vector<t_real>;
using t_vec_cplx = ublas::vector<t_cplx>;

#include "core/skx_default_gs.cxx"


const t_real g_eps = 1e-5;
const t_real g_weight_eps = 1e-6;

#define E_BINS 200
#define NUM_ANGLES 512
#define COL_SIZE 15


void calc_disp(const t_vec& Gvec,
	const t_vec& Bvec, const t_vec& Pvec,
	t_real q, t_real q_oop, bool bProj = true)
{
	// vector perpendicular to the pinning, but in the skyrmion plane
	t_vec Pperpvec = tl2::cross_3(Pvec, Bvec);
	Pperpvec /= tl2::veclen(Pperpvec);

	Skx<t_real, t_cplx, DEF_SKX_ORDER> skx;
	Heli<t_real, t_cplx, DEF_HELI_ORDER> heli;

	skx.SetFourier(_get_skx_gs<t_vec_cplx>());

	skx.SetProjNeutron(bProj);
	heli.SetProjNeutron(bProj);

	t_real T = 28.5;
	skx.SetT(-1000., false);
	heli.SetT(-1000., false);
	heli.SetT(T, true);

	t_real bc2 = skx.GetBC2(false);
	skx.SetB(bc2/2., false);
	heli.SetB(bc2/2., false);
	//heli.SetB(0.17, true);
	heli.SetB(0.158, true);

	skx.SetFilterZeroWeight(1);
	heli.SetFilterZeroWeight(1);

	skx.SetWeightEps(g_weight_eps);
	heli.SetWeightEps(g_weight_eps);

	skx.SetCoords(Bvec[0], Bvec[1], Bvec[2], Pvec[0], Pvec[1], Pvec[2]);
	heli.SetCoords(Bvec[0], Bvec[1], Bvec[2]);

	skx.SetG(Gvec[0], Gvec[1], Gvec[2]);
	heli.SetG(Gvec[0], Gvec[1], Gvec[2]);

	t_real Erange = 0.1;

	t_real angle_begin = -45/180.*M_PI;
	t_real angle_end = 270/180.*M_PI - 45/180.*M_PI;
	t_real angle_delta = 2*M_PI/t_real(NUM_ANGLES);

	auto histWeightsNSF = hist::make_histogram(hist::axis::regular<t_real>(E_BINS, -Erange, Erange, "E"));
	auto histWeightsSF = hist::make_histogram(hist::axis::regular<t_real>(E_BINS, -Erange, Erange, "E"));

	auto histWeightsHeliNSF = hist::make_histogram(hist::axis::regular<t_real>(E_BINS, -Erange, Erange, "E"));
	auto histWeightsHeliSF = hist::make_histogram(hist::axis::regular<t_real>(E_BINS, -Erange, Erange, "E"));


	std::ofstream ofstr_raw("weightsum_skx.dat");
	std::ofstream ofstr_raw_heli("weightsum_heli.dat");

	// write file header
	for(std::ostream* ostr : {&ofstr_raw, &ofstr_raw_heli})
	{
		ostr->precision(8);
		(*ostr)  << std::left << std::setw(COL_SIZE) << "# angle"
			<< " " << std::left << std::setw(COL_SIZE) << "qh"
			<< " " << std::left << std::setw(COL_SIZE) << "qk"
			<< " " << std::left << std::setw(COL_SIZE) << "ql"
			<< " " << std::left << std::setw(COL_SIZE) << "E"
			<< " " << std::left << std::setw(COL_SIZE) << "wSF1"
			<< " " << std::left << std::setw(COL_SIZE) << "wSF2"
			<< " " << std::left << std::setw(COL_SIZE) << "wNSF"
			<< "\n";
	}


	for(t_real angle=angle_begin; angle<angle_end; angle+=angle_delta)
	{
		t_vec qvec = q * (Pvec*std::cos(angle) + Pperpvec*std::sin(angle));  // in skx plane
		qvec += q_oop * Bvec;                                                // out of skx plane
		t_vec Qvec = Gvec + qvec;

		std::cout << "# angle: " << angle/M_PI*180.
			<< ", Q = (" << Qvec[0] << ", " << Qvec[1] << ", " << Qvec[2] << ")"
			<< ", |q| = " << tl2::veclen(qvec) << "."
			<< std::endl;

		{
			// skyrmion dispersion
			auto [Es, wsUnpol, wsSF1, wsSF2, wsNSF] = skx.GetDisp(Qvec[0], Qvec[1], Qvec[2], -Erange, Erange);
			for(std::size_t i=0; i<Es.size(); ++i)
			{
				histWeightsNSF(Es[i], hist::weight(wsNSF[i]*0.5));
				histWeightsSF(Es[i], hist::weight(wsSF1[i]));

				ofstr_raw << std::left << std::setw(COL_SIZE) << angle
					<< " " << std::left << std::setw(COL_SIZE) << qvec[0]
					<< " " << std::left << std::setw(COL_SIZE) << qvec[1]
					<< " " << std::left << std::setw(COL_SIZE) << qvec[2]
					<< " " << std::left << std::setw(COL_SIZE) << Es[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsSF1[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsSF2[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsNSF[i]
					<< std::endl;
			}
		}

		{
			// conical dispersion
			auto [EsH, wsUnpolH, wsSF1H, wsSF2H, wsNSFH] = heli.GetDisp(Qvec[0], Qvec[1], Qvec[2], -Erange, Erange);
			for(std::size_t i=0; i<EsH.size(); ++i)
			{
				histWeightsHeliNSF(EsH[i], hist::weight(wsNSFH[i]*0.5));
				histWeightsHeliSF(EsH[i], hist::weight(wsSF1H[i]));

				ofstr_raw_heli << std::left << std::setw(COL_SIZE) << angle
					<< " " << std::left << std::setw(COL_SIZE) << qvec[0]
					<< " " << std::left << std::setw(COL_SIZE) << qvec[1]
					<< " " << std::left << std::setw(COL_SIZE) << qvec[2]
					<< " " << std::left << std::setw(COL_SIZE) << EsH[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsSF1H[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsSF2H[i]
					<< " " << std::left << std::setw(COL_SIZE) << wsNSFH[i]
					<< std::endl;
			}
		}
	}


	std::ofstream ofstrBinned("weightbin.dat");
	ofstrBinned.precision(8);

	ofstrBinned << std::left << std::setw(COL_SIZE) << "# E (meV)"
		<< " " << std::left << std::setw(COL_SIZE) << "bose"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_skx_sf"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_skx_nsf"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_skx_sum"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_skx_sum_bose"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_heli_sf"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_heli_nsf"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_heli_sum"
		<< " " << std::left << std::setw(COL_SIZE) << "weight_heli_sum_bose"
		<< "\n";

	auto iterHeliNSF = boost::histogram::indexed(histWeightsHeliNSF).begin();
	auto iterHeliSF = boost::histogram::indexed(histWeightsHeliSF).begin();
	auto iterNSF = boost::histogram::indexed(histWeightsNSF).begin();
	for(const auto& val : boost::histogram::indexed(histWeightsSF))
	{
		t_real E = val.bin().lower() + 0.5*(val.bin().upper() - val.bin().lower());
		t_real w = *val / t_real{E_BINS};

		t_real E_nsf = iterNSF->bin().lower() + 0.5*(iterNSF->bin().upper() - iterNSF->bin().lower());
		t_real w_nsf = **iterNSF / t_real{E_BINS};

		t_real E_h = iterHeliSF->bin().lower() + 0.5*(iterHeliSF->bin().upper() - iterHeliSF->bin().lower());
		t_real w_h = **iterHeliSF / t_real{E_BINS};

		t_real E_h_nsf = iterHeliNSF->bin().lower() + 0.5*(iterHeliNSF->bin().upper() - iterHeliNSF->bin().lower());
		t_real w_h_nsf = **iterHeliNSF / t_real{E_BINS};

		if(!tl2::float_equal<t_real>(E, E_h, g_eps) || !tl2::float_equal<t_real>(E, E_h_nsf, g_eps) ||
			!tl2::float_equal<t_real>(E, E_nsf, g_eps))

		{
			std::cerr << "Energy binning mismatch: " << E << " != " << E_h << std::endl;
			break;
		}

		t_real bose = tl2::bose(E, T);

		ofstrBinned << std::left << std::setw(COL_SIZE) << E
			<< " " << std::left << std::setw(COL_SIZE) << bose
			<< " " << std::left << std::setw(COL_SIZE) << w
			<< " " << std::left << std::setw(COL_SIZE) << w_nsf
			<< " " << std::left << std::setw(COL_SIZE) << (w + w_nsf)
			<< " " << std::left << std::setw(COL_SIZE) << (w + w_nsf)*bose
			<< " " << std::left << std::setw(COL_SIZE) << w_h
			<< " " << std::left << std::setw(COL_SIZE) << w_h_nsf
			<< " " << std::left << std::setw(COL_SIZE) << (w_h + w_h_nsf)
			<< " " << std::left << std::setw(COL_SIZE) << (w_h + w_h_nsf)*bose
			<< std::endl;

		++iterHeliSF;
		++iterHeliNSF;
		++iterNSF;
	}
}


int main()
{
	t_vec Gvec = tl2::make_vec<t_vec>({ 0., 0., 0. });  // lattice vector
	t_vec Pvec = tl2::make_vec<t_vec>({ 1., 1., 0. });  // pinning vector
	t_vec Bvec = tl2::make_vec<t_vec>({ 0., 0., 1. });  // field direction

	Pvec /= tl2::veclen(Pvec);
	Bvec /= tl2::veclen(Bvec);

	t_real q = 0.0123;    // momentum transfer in skyrmion plane
	t_real q_oop = 0.0;   // momentum transfer out of skyrmion plane
	bool proj = true;     // using the orthogonal 1-|Q><Q| projector

	calc_disp(Gvec, Bvec, Pvec, q, q_oop, proj);
	return 0;
}
