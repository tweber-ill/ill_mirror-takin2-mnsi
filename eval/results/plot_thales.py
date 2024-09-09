#!/usr/bin/python3
#
# plots thales data files
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
def load_data(datfile, mergefiles = [], I_scale = 1.):
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
		Is.append(I * I_scale)
		Is_err.append(I_err * I_scale)

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
outdir = "./thales/"

try:
	os.makedirs(outdir)
except FileExistsError:
	pass

import matplotlib.pyplot as mplt
mplt.rcParams["font.size"] = 16

ms = 6
cs = 6
inc_sig_14 = 0.117818 / 2. / (2. * numpy.log(2.))**0.5
inc_amp_14 = 1.2

inc_sig_12 = 0.0626775 / 2. / (2. * numpy.log(2.))**0.5
inc_amp_12 = 0.55

I_scale = 1e4
def S_scale(T):
	#return I_scale
	return (-0.0625 * T + 2.5) * I_scale


# -----------------------------------------------------------------------------
print("Generating plot 1...")
fig, (plt) = mplt.subplots(1, 1)

plt.set_xlim(-0.175, 1.075)
plt.set_ylim(0, 1.)

plt.set_xlabel("E (meV)")
plt.set_ylabel("S (a.u.)")

plt.annotate("scan (v), Thales\n(1.07 0.93 0)", xy=(0.62, 0.53), xycoords="axes fraction")

scale = 2.0  # 1.4
offs = 5e-6

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/ill_thales/exp_INTER-477/rawdata/024813", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="o", markersize=ms, capsize=cs, ls="none", color=get_col(0,3), label="%.1f K (coni)" % T)
if os.path.exists("coni_107_perp_Bvert_25K.dat"):
	convo = numpy.loadtxt("coni_107_perp_Bvert_25K.dat")
	#plt.plot(convo[:,3], convo[:,5]*S_scale(T) + gauss(convo[:,3], 0., inc_sig_14, inc_amp_14, 0.), color=get_col(0,3))
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_14, inc_amp_14, 0.), color=get_col(0,3))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/ill_thales/exp_INTER-477/rawdata/024808", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="s", markersize=ms, capsize=cs, ls="none", color=get_col(1,3), label="%.1f K (coni)" % T)
if os.path.exists("coni_107_perp_Bvert_27K.dat"):
	convo = numpy.loadtxt("coni_107_perp_Bvert_27K.dat")
	#plt.plot(convo[:,3], convo[:,5]*S_scale(T) + gauss(convo[:,3], 0., inc_sig_14, inc_amp_14, 0.), color=get_col(1,3))
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_14, inc_amp_14, 0.), color=get_col(1,3))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/ill_thales/exp_INTER-477/rawdata/024801", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="^", markersize=ms, capsize=cs, ls="none", color=get_col(2,3), label="%.1f K (skx)" % T)
if os.path.exists("skx_107_perp_Bvert.dat"):
	convo = numpy.loadtxt("skx_107_perp_Bvert.dat")
	#plt.plot(convo[:,3], convo[:,5]*S_scale(T) + gauss(convo[:,3], 0., inc_sig_14, inc_amp_14, 0.), color=get_col(2,3))
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_14, inc_amp_14, 0.), color=get_col(2,3))

plt.legend()
fig.tight_layout()
fig.savefig(outdir + "107_perp_Bvert.pdf")
print("Saved plot 1 as 107_perp_Bvert.pdf.\n")
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
print("Generating plot 2...")
fig, (plt) = mplt.subplots(1, 1)

plt.set_xlim(-0.175, 1.075)
plt.set_ylim(0, 0.7)

plt.set_xlabel("E (meV)")
plt.set_ylabel("S (a.u.)")

plt.annotate("scan (iv), Thales\n(1.06 0.94 0)", xy=(0.62, 0.5), xycoords="axes fraction")

scale = 1.15
offs = 5e-6
num_scans = 4
idx = 0
show_25K = False
if not show_25K:
	num_scans = 3

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/ill_thales/exp_INTER-477/rawdata/024828", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="o", markersize=ms, capsize=cs, ls="none", color=get_col(idx, num_scans), label="%.1f K (coni)" % T)
if os.path.exists("coni_106_perp_Bvert_20K.dat"):
	convo = numpy.loadtxt("coni_106_perp_Bvert_20K.dat")
	#plt.plot(convo[:,3], convo[:,5]*S_scale(T) + gauss(convo[:,3], 0., inc_sig_12, inc_amp_12, 0.), color=get_col(0,3))
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_12, inc_amp_12, 0.), color=get_col(idx, num_scans))
idx += 1

if show_25K:
	(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/ill_thales/exp_INTER-477/rawdata/024823", I_scale=I_scale,
		mergefiles=["../data/ill_thales/exp_INTER-477/rawdata/024793"])
	plt.errorbar(Es, Is, Is_err, marker="s", markersize=ms, capsize=cs, ls="none", color=get_col(idx, num_scans), label="%.1f K (coni)" % T)
	if os.path.exists("coni_106_perp_Bvert_25K.dat"):
		convo = numpy.loadtxt("coni_106_perp_Bvert_25K.dat")
		#plt.plot(convo[:,3], convo[:,5]*S_scale(T) + gauss(convo[:,3], 0., inc_sig_12, inc_amp_12, 0.), color=get_col(1,3))
		plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_12, inc_amp_12, 0.), color=get_col(idx, num_scans))
	idx += 1

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/ill_thales/exp_INTER-477/rawdata/024788", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="^", markersize=ms, capsize=cs, ls="none", color=get_col(idx, num_scans), label="%.1f K (coni)" % T)
if os.path.exists("coni_106_perp_Bvert_27K.dat"):
	convo = numpy.loadtxt("coni_106_perp_Bvert_27K.dat")
	#plt.plot(convo[:,3], convo[:,5]*S_scale(T) + gauss(convo[:,3], 0., inc_sig_12, inc_amp_12, 0.), color=get_col(0,3))
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_12, inc_amp_12, 0.), color=get_col(idx, num_scans))
idx += 1

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/ill_thales/exp_INTER-477/rawdata/024816", I_scale=I_scale,
	mergefiles=["../data/ill_thales/exp_INTER-477/rawdata/024781"])
plt.errorbar(Es, Is, Is_err, marker="v", markersize=ms, capsize=cs, ls="none", color=get_col(idx, num_scans), label="%.1f K (skx)" % T)
if os.path.exists("skx_106_perp_Bvert.dat"):
	convo = numpy.loadtxt("skx_106_perp_Bvert.dat")
	#plt.plot(convo[:,3], convo[:,5]*S_scale(T) + gauss(convo[:,3], 0., inc_sig_12, inc_amp_12, 0.), color=get_col(2,3))
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_12, inc_amp_12, 0.), color=get_col(idx, num_scans))

plt.legend()
fig.tight_layout()
fig.savefig(outdir + "106_perp_Bvert.pdf")
print("Saved plot 2 as 106_perp_Bvert.pdf.\n")
# -----------------------------------------------------------------------------
