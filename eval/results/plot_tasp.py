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

I_scale = 1e3
def S_scale(T):
	#return I_scale * 0.9
	return (-0.0625 * T + 2.5) * I_scale


# -----------------------------------------------------------------------------
print("Generating plot 1...")
fig, (plt) = mplt.subplots(1, 1)

plt.set_xlim(-0.625, 0.625)
plt.set_ylim(0, 4)

plt.set_xlabel("E (meV)")
plt.set_ylabel("S (a.u.)")

plt.annotate("scan (i), TASP\n(0.944 0.944 0)", xy=(0.6, 0.85), xycoords="axes fraction")
plt.annotate("B = -195 mT", xy=(0.1, 0.5), xycoords="axes fraction")

# arrows
# ./weight --dyntype=h --Gx=1 --Gy=1 --Gz=0 --Bx=-1 --By=-1 --Bz=0 --Px=1 --Py=-1 --Pz=0 --Qx=0.944 --Qy=0.944 --Qz=0 --T=18.03 --B=0.169845
# ./weight --dyntype=h --Gx=1 --Gy=1 --Gz=0 --Bx=-1 --By=-1 --Bz=0 --Px=1 --Py=-1 --Pz=0 --Qx=0.944 --Qy=0.944 --Qz=0 --T=24.04 --B=0.169845
# ./weight --dyntype=h --Gx=1 --Gy=1 --Gz=0 --Bx=-1 --By=-1 --Bz=0 --Px=1 --Py=-1 --Pz=0 --Qx=0.944 --Qy=0.944 --Qz=0 --T=27.01 --B=0.169845
# ./weight --dyntype=s --Gx=1 --Gy=1 --Gz=0 --Bx=-1 --By=-1 --Bz=0 --Px=1 --Py=-1 --Pz=0 --Qx=0.944 --Qy=0.944 --Qz=0 --T=28.52 --B=0.169845
arrow_len = 0.35
plt.annotate("18K (c)", xytext=(0.289, 1.4), xy=(0.289, 1.4-arrow_len), color=get_col(0,4), \
	arrowprops={"color" : get_col(0,4), "width":2., "headwidth":8., "headlength":8., "shrink":0.1})
plt.annotate("24K (c)", xytext=(0.221, 1.9), xy=(0.221, 1.9-arrow_len), color=get_col(1,4), \
	arrowprops={"color" : get_col(1,4), "width":2., "headwidth":8., "headlength":8., "shrink":0.1})
plt.annotate("27K (c)", xytext=(0.174, 2.5), xy=(0.174, 2.5-arrow_len), color=get_col(2,4), \
	arrowprops={"color" : get_col(2,4), "width":2., "headwidth":8., "headlength":8., "shrink":0.1})
plt.annotate("28.5K (s)", xytext=(0.132, 3.1), xy=(0.132, 3.1-arrow_len), color=get_col(3,4), \
	arrowprops={"color" : get_col(3,4), "width":2., "headwidth":8., "headlength":8., "shrink":0.1})

scale = 4
offs = 1e-4
(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003343.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="o", markersize=ms, capsize=cs, ls="none", color=get_col(0,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_18K.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_18K.dat")
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(0,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003219.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="s", markersize=ms, capsize=cs, ls="none", color=get_col(1,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_24K.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_24K.dat")
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(1,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003225.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="^", markersize=ms, capsize=cs, ls="none", color=get_col(2,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_27K.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_27K.dat")
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(2,4))

#scale = 3.5
(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003238.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="v", markersize=ms, capsize=cs, ls="none", color=get_col(3,4), label="%.1f K (skx)" % T)
if os.path.exists("convo_tasp_0944_para.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para.dat")
	#plt.plot(convo[:,3], convo[:,5]*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(3,4))
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(3,4))

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

plt.annotate("scan (i), TASP\n(0.944 0.944 0)", xy=(0.6, 0.85), xycoords="axes fraction")

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003238.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="v", markersize=ms, capsize=cs, ls="none", color=get_col(0,4,2), label="%.1f K (skx)" % T)
if os.path.exists("convo_tasp_0944_para.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para.dat")
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(3,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003242.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="^", markersize=ms, capsize=cs, ls="none", color=get_col(1,4,2), label="%.1f K (fd)" % T)

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003251.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="s", markersize=ms, capsize=cs, ls="none", color=get_col(2,4,2), label="%.1f K (fd)" % T)

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003255.dat", I_scale=I_scale)
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
plt.set_ylim(0, 4)

plt.set_xlabel("E (meV)")
plt.set_ylabel("S (a.u.)")

plt.annotate("scan (i), TASP\n(0.944 0.944 0)", xy=(0.05, 0.85), xycoords="axes fraction")
plt.annotate("B = 195 mT", xy=(0.05, 0.7), xycoords="axes fraction")

# arrows
# ./weight --dyntype=h --Gx=1 --Gy=1 --Gz=0 --Bx=1 --By=1 --Bz=0 --Px=1 --Py=-1 --Pz=0 --Qx=0.944 --Qy=0.944 --Qz=0 --T=18.03 --B=0.169845
# ./weight --dyntype=h --Gx=1 --Gy=1 --Gz=0 --Bx=1 --By=1 --Bz=0 --Px=1 --Py=-1 --Pz=0 --Qx=0.944 --Qy=0.944 --Qz=0 --T=24.04 --B=0.169845
# ./weight --dyntype=h --Gx=1 --Gy=1 --Gz=0 --Bx=1 --By=1 --Bz=0 --Px=1 --Py=-1 --Pz=0 --Qx=0.944 --Qy=0.944 --Qz=0 --T=27.01 --B=0.169845
# ./weight --dyntype=s --Gx=1 --Gy=1 --Gz=0 --Bx=1 --By=1 --Bz=0 --Px=1 --Py=-1 --Pz=0 --Qx=0.944 --Qy=0.944 --Qz=0 --T=28.52 --B=0.169845
arrow_len = 0.35
plt.annotate("", xytext=(-0.289, 1.2), xy=(-0.289, 1.2-arrow_len), color=get_col(0,4), \
	arrowprops={"color" : get_col(0,4), "width":2., "headwidth":8., "headlength":8., "shrink":0.1})
plt.annotate("", xytext=(-0.221, 1.6), xy=(-0.221, 1.6-arrow_len), color=get_col(1,4), \
	arrowprops={"color" : get_col(1,4), "width":2., "headwidth":8., "headlength":8., "shrink":0.1})
plt.annotate("", xytext=(-0.174, 2.1), xy=(-0.174, 2.1-arrow_len), color=get_col(2,4), \
	arrowprops={"color" : get_col(2,4), "width":2., "headwidth":8., "headlength":8., "shrink":0.1})
plt.annotate("", xytext=(-0.132, 2.5), xy=(-0.132, 2.5-arrow_len), color=get_col(3,4), \
	arrowprops={"color" : get_col(3,4), "width":2., "headwidth":8., "headlength":8., "shrink":0.1})

label_offs_x = -0.225
label_offs_y = -0.2
plt.annotate("18K (c)", xy=(-0.289+label_offs_x, 1.2+label_offs_y), color=get_col(0,4))
plt.annotate("24K (c)", xy=(-0.221+label_offs_x, 1.6+label_offs_y), color=get_col(1,4))
plt.annotate("27K (c)", xy=(-0.174+label_offs_x, 2.1+label_offs_y), color=get_col(2,4))
plt.annotate("28.5K (s)", xy=(-0.132+label_offs_x-0.05, 2.5+label_offs_y), color=get_col(3,4))

scale = 4
offs = 1e-4
(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003339.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="o", markersize=ms, capsize=cs, ls="none", color=get_col(0,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_18K_Bflipped.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_18K_Bflipped.dat")
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(0,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003222.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="s", markersize=ms, capsize=cs, ls="none", color=get_col(1,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_24K_Bflipped.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_24K_Bflipped.dat")
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(1,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003228.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="^", markersize=ms, capsize=cs, ls="none", color=get_col(2,4), label="%.1f K (coni)" % T)
if os.path.exists("convo_tasp_0944_para_27K_Bflipped.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_27K_Bflipped.dat")
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(2,4))

#scale = 3.5
(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003233.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="v", markersize=ms, capsize=cs, ls="none", color=get_col(3,4), label="%.1f K (skx)" % T)
if os.path.exists("convo_tasp_0944_para_Bflipped.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_Bflipped.dat")
	plt.plot(convo[:,3], (convo[:,4]*scale + offs)*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(3,4))

plt.legend(handletextpad=0.)
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

plt.annotate("scan (i), TASP\n(0.944 0.944 0)", xy=(0.6, 0.85), xycoords="axes fraction")

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003233.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="v", markersize=ms, capsize=cs, ls="none", color=get_col(0,4,2), label="%.1f K (skx)" % T)
if os.path.exists("convo_tasp_0944_para_Bflipped.dat"):
	convo = numpy.loadtxt("convo_tasp_0944_para_Bflipped.dat")
	plt.plot(convo[:,3], convo[:,5]*S_scale(T) + gauss(convo[:,3], 0., inc_sig_145, inc_amp_145, 0.), color=get_col(3,4))

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003245.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="^", markersize=ms, capsize=cs, ls="none", color=get_col(1,4,2), label="%.1f K (fd)" % T)

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003248.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="s", markersize=ms, capsize=cs, ls="none", color=get_col(2,4,2), label="%.1f K (fd)" % T)

(hs, ks, ls, Es, Is, Is_err, T) = load_data("../data/psi_tasp/exp_20181324_2/tasp2018n003258.dat", I_scale=I_scale)
plt.errorbar(Es, Is, Is_err, marker="o", markersize=ms, capsize=cs, ls="none", color=get_col(3,4,2), label="%.1f K (pm)" % T)

plt.legend()
fig.tight_layout()
fig.savefig(outdir + "0944_para_pB_above_skx.pdf")
print("Saved plot 4 as 0944_para_pB_above_skx.pdf.\n")
# -----------------------------------------------------------------------------
