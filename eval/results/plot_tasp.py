#!/usr/bin/python3
#
# plots tasp data files
# @author Tobias Weber <tweber@ill.fr>
# @date apr-2024
# @license see 'LICENSE' file
#

import os
import sys
import numpy
import scipy.constants as co
import instr


#
# normalise detector (y) by monitor (m) counts
# y_new = y / m
# dy_new = 1/m dy - y/m^2 dm
#
def norm_counts_to_mon(y, dy, m, dm):
	val = y / m;
	err = ((dy/m)**2. + (dm*y/(m*m))**2.)**0.5
	return [val, err]


#
# get energy transfer from ki and kf
#
def get_E(ki, kf):
	E_to_k2 = 2.*co.neutron_mass/co.hbar**2. * co.elementary_charge/1000. * 1e-20
	return (ki**2. - kf**2.) / E_to_k2


#
# load an instrument data file
#
def load_data(datfile, mergefiles = [], Iscale = 1.):
	print("Loading \"%s\"." % (datfile))
	dat = instr.FileInstrBaseD.LoadInstr(datfile)
	if dat == None:
		return

	for mergefile in mergefiles:
		tomerge = instr.FileInstrBaseD.LoadInstr(mergefile)
		if tomerge != None:
			print("Merging with \"%s\"." % (mergefile))
			dat.MergeWith(tomerge)


	hs = []; ks = []; ls = []; Es = []
	Is = []; Is_err = []

	cntcol = dat.GetCol(dat.GetCountVar())
	moncol = dat.GetCol(dat.GetMonVar())
	T = numpy.mean(dat.GetCol("TT"))

	# iterate scan points
	for point_idx in range(cntcol.size()):
		(h, k, l, ki, kf) = dat.GetScanHKLKiKf(point_idx)
		E = get_E(ki, kf)

		hs.append(h)
		ks.append(k)
		ls.append(l)
		Es.append(E)

		counts = cntcol[point_idx]
		mon = moncol[point_idx]

		if counts == 0:
			counts_err = 1.
		else:
			counts_err = counts**0.5
		if mon == 0:
			mon_err = 1.
		else:
			mon_err = mon**0.5

		[I, I_err] = norm_counts_to_mon(counts, counts_err, mon, mon_err)
		Is.append(I * Iscale)
		Is_err.append(I_err * Iscale)

	return (hs, ks, ls, Es, Is, Is_err, T)


#
# gaussian
#
def gauss(x, x0, sigma, amp, offs):
	norm = 1. / ((2.*numpy.pi)**0.5 * sigma)
	return amp * norm * numpy.exp(-0.5 * ((x - x0)/sigma)**2.) + offs


#
# lorentzian
#
def lorentz(x, x0, hwhm, amp, offs):
	return amp * hwhm**2. / ((x - x0)**2. + hwhm**2.) + offs


#
# get colour string
#
def get_col(idx, num_cols, style = 1):
	if style == 1:
		r = int(255. * float(idx) / float(num_cols - 1))
		g = 0
		b = 255 - r
	elif style == 2:
		r = 255 - int(255. * float(idx) / float(num_cols - 1))
		g = 0
		b = 0

	return "#%02x%02x%02x" % (r, g, b)


# =============================================================================
outdir = "./tasp/"

try:
	os.makedirs(outdir)
except FileExistsError:
	pass

import matplotlib.pyplot as mplt
mplt.rcParams["font.size"] = 16

ms = 6
cs = 6
inc_sig_145 = 0.163753 / 2. / (2. * numpy.log(2.))**0.5
inc_amp_145 = 0.5

Iscale = 1e3
def Iscale_func(T):
	return Iscale * 0.9
	#return (-0.0625 * T + 2.5) * Iscale


# -----------------------------------------------------------------------------
print("Generating plot 1...")
fig, (plt) = mplt.subplots(1, 1)

plt.set_xlim(-0.625, 0.625)
#plt.set_ylim(0, 2)

plt.set_xlabel("E (meV)")
plt.set_ylabel("S (a.u.)")

plt.annotate("(0.944 0.944 0)", xy=(0.6, 0.85), xycoords="axes fraction")

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003343.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="o", markersize=ms, capsize=cs, ls="none", color=get_col(0,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_18K.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_18K.dat")
	plt.plot(convo[:,3], convo[:,5]*Iscale_func(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(0,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003219.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="s", markersize=ms, capsize=cs, ls="none", color=get_col(1,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_24K.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_24K.dat")
	plt.plot(convo[:,3], convo[:,5]*Iscale_func(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(1,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003225.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="^", markersize=ms, capsize=cs, ls="none", color=get_col(2,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_27K.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_27K.dat")
	plt.plot(convo[:,3], convo[:,5]*Iscale_func(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(2,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003238.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="v", markersize=ms, capsize=cs, ls="none", color=get_col(3,4), label="%.1f K (skx)" % T)
if os.path.exists("convo_tasp_0944_para.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para.dat")
	plt.plot(convo[:,3], convo[:,5]*Iscale_func(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(3,4))

plt.legend()
fig.tight_layout()
fig.savefig(outdir + "0944_para_mB.pdf")
print("Saved plot 1 as 0944_para_mB.pdf.\n")
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
print("Generating plot 2...")
fig, (plt) = mplt.subplots(1, 1)

plt.set_xlim(-0.625, 0.625)
#plt.set_ylim(0, 2)

plt.set_xlabel("E (meV)")
plt.set_ylabel("S (a.u.)")

plt.annotate("(0.944 0.944 0)", xy=(0.6, 0.85), xycoords="axes fraction")

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003238.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="v", markersize=ms, capsize=cs, ls="none", color=get_col(0,4,2), label="%.1f K (skx)" % T)
if os.path.exists("convo_tasp_0944_para.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para.dat")
	plt.plot(convo[:,3], convo[:,5]**Iscale_func(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(3,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003242.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="^", markersize=ms, capsize=cs, ls="none", color=get_col(1,4,2), label="%.1f K (fd)" % T)

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003251.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="s", markersize=ms, capsize=cs, ls="none", color=get_col(2,4,2), label="%.1f K (fd)" % T)

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003255.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="o", markersize=ms, capsize=cs, ls="none", color=get_col(3,4,2), label="%.1f K (pm)" % T)

plt.legend()
fig.tight_layout()
fig.savefig(outdir + "0944_para_mB_above_skx.pdf")
print("Saved plot 2 as 0944_para_mB_above_skx.pdf.\n")
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
print("Generating plot 3...")
fig, (plt) = mplt.subplots(1, 1)

plt.set_xlim(-0.625, 0.625)
#plt.set_ylim(0, 2)

plt.set_xlabel("E (meV)")
plt.set_ylabel("S (a.u.)")

plt.annotate("(0.944 0.944 0)", xy=(0.6, 0.85), xycoords="axes fraction")

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003339.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="o", markersize=ms, capsize=cs, ls="none", color=get_col(0,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_18K_Bflipped.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_18K_Bflipped.dat")
	plt.plot(convo[:,3], convo[:,5]**Iscale_func(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(0,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003222.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="s", markersize=ms, capsize=cs, ls="none", color=get_col(1,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_24K_Bflipped.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_24K_Bflipped.dat")
	plt.plot(convo[:,3], convo[:,5]**Iscale_func(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(1,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003228.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="^", markersize=ms, capsize=cs, ls="none", color=get_col(2,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_27K_Bflipped.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_27K_Bflipped.dat")
	plt.plot(convo[:,3], convo[:,5]**Iscale_func(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(2,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003233.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="v", markersize=ms, capsize=cs, ls="none", color=get_col(3,4), label="%.1f K (skx)" % T)
if os.path.exists("convo_tasp_0944_para_Bflipped.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_Bflipped.dat")
	plt.plot(convo[:,3], convo[:,5]**Iscale_func(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(3,4))

plt.legend()
fig.tight_layout()
fig.savefig(outdir + "0944_para_pB.pdf")
print("Saved plot 3 as 0944_para_pB.pdf.\n")
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
print("Generating plot 4...")
fig, (plt) = mplt.subplots(1, 1)

plt.set_xlim(-0.625, 0.625)
#plt.set_ylim(0, 2)

plt.set_xlabel("E (meV)")
plt.set_ylabel("S (a.u.)")

plt.annotate("(0.944 0.944 0)", xy=(0.6, 0.85), xycoords="axes fraction")

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003233.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="v", markersize=ms, capsize=cs, ls="none", color=get_col(0,4,2), label="%.1f K (skx)" % T)
if os.path.exists("convo_tasp_0944_para_Bflipped.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_Bflipped.dat")
	plt.plot(convo[:,3], convo[:,5]**Iscale_func(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(3,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003245.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="^", markersize=ms, capsize=cs, ls="none", color=get_col(1,4,2), label="%.1f K (fd)" % T)

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003248.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="s", markersize=ms, capsize=cs, ls="none", color=get_col(2,4,2), label="%.1f K (fd)" % T)

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003258.dat", Iscale=Iscale)
plt.errorbar(Es, Is, Is_err, marker="o", markersize=ms, capsize=cs, ls="none", color=get_col(3,4,2), label="%.1f K (pm)" % T)

plt.legend()
fig.tight_layout()
fig.savefig(outdir + "0944_para_pB_above_skx.pdf")
print("Saved plot 4 as 0944_para_pB_above_skx.pdf.\n")
# -----------------------------------------------------------------------------