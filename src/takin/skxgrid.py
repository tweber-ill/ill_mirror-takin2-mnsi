#
# skx grid
# @author tweber@ill.fr
# @date 4-feb-2020
# @license GPLv2 (see 'LICENSE' file)
#

import numpy as np


# -----------------------------------------------------------------------------
# Test files and configuration

use_thales = True
use_reseda = False
slice = 1

if use_thales:
    idxfile = "/home/tw/tmp/skx/skxdyn_asym.idx"
    datafile = "/home/tw/tmp/skx/skxdyn_asym.bin"

    hstep = 0.001
    hmin = -0.096
    hmax = 0.096

    kstep = 0.001
    kmin = -0.096
    kmax = 0.096

    lstep = 0.001
    lmin = -0.096
    lmax = 0.096

    xlim1 = -0.09
    xlim2 = 0.09
    ylim1 = -1.
    ylim2 = 1.


if use_reseda:
    idxfile = "/home/tw/tmp/skx/skxdyn_reseda.idx"
    datafile = "/home/tw/tmp/skx/skxdyn_reseda.bin"

    hstep = 0.0006
    hmin = -0.03
    hmax = 0.03 + hstep

    kstep = 0.0006
    kmin = -0.03
    kmax = 0.03 + kstep

    lstep = 0.0006
    lmin = -0.03
    lmax = 0.03 + lstep

    xlim1 = -0.03
    xlim2 = 0.03
    ylim1 = -0.3
    ylim2 = 0.3


numvals = 4     # [E, w1, w2, w3]
plot_hklsteps = 512


# selection of slices
if slice==0:
    # longitudinal
    plot_hklbegin = np.array([-0.09, -0.09, -0.])
    plot_hklend = np.array([0.09, 0.09, 0.])
    plot_dir = 0
elif slice==1:
    # transversal
    plot_hklbegin = np.array([-0.09, 0.09, 0.])
    plot_hklend = np.array([0.09, -0.09, 0.])
    plot_dir = 0
elif slice==2:
    # up
    plot_hklbegin = np.array([0., 0., -0.09])
    plot_hklend = np.array([0., 0., 0.09])
    plot_dir = 2
elif slice==3:
    plot_hklbegin = np.array([0., -0.03, 0.])
    plot_hklend = np.array([0., 0.03, 0.])
    plot_dir = 1
elif slice==4:
    plot_hklbegin = np.array([0.03, -0.03, -0.])
    plot_hklend = np.array([-0.03, 0.03, 0.])
    plot_dir = 0
# -----------------------------------------------------------------------------



# using q, not Q
def hkl_to_idx(hkl):
    [h, k, l] = hkl

    # clamp values to boundaries
    if h < hmin:
        h = hmin
    if k < kmin:
        k = kmin
    if l < lmin:
        l = lmin
    if h >= hmax:
        h = hmax - hstep
    if k >= kmax:
        k = kmax - kstep
    if l >= lmax:
        l = lmax - lstep

    # max dimensions
    iHSize = int(round((hmax - hmin) / hstep))
    iKSize = int(round((kmax - kmin) / kstep))
    iLSize = int(round((lmax - lmin) / lstep))

    # position indices
    iH = int(round(((h - hmin) / hstep)))
    iK = int(round(((k - kmin) / kstep)))
    iL = int(round(((l - lmin) / lstep)))

    # clamp again
    if iH >= iHSize:
        iH = iHSize-1
    if iK >= iKSize:
        iK = iKSize-1
    if iL >= iLSize:
        iL = iLSize-1
    if iH < 0:
        iH = 0
    if iK < 0:
        iK = 0
    if iL < 0:
        iL = 0

    return int(iH*iKSize*iLSize + iK*iLSize + iL)



def getE(idxfilehandle, datafilehandle, hkl):
    hklidx = hkl_to_idx(hkl)
    #print("Index into index file: %d." % hklidx)

    idx = np.memmap(idxfilehandle, dtype="uint64", mode="r", offset=int(hklidx*8))[0]
    #print("Index into data file: %d." % idx)

    data = np.memmap(datafilehandle, dtype="uint32", mode="r", offset=int(idx))
    numbranches = data[0]
    #print("Number of dispersion branches: %d." % numbranches)

    values = np.ndarray.view(data[1 : 1+numbranches*numvals*2], dtype="float64")
    values.shape = (numbranches, numvals)
    return values




# plot an example dispersion
def plot_disp(idxfilehandle, datafilehandle, hklbegin, hklend, hklsteps):
    import matplotlib.pyplot as plot
    symscale = 2.
    eps = 1e-8

    qs_h = []
    qs_k = []
    qs_l = []
    Es = []
    wsSF1 = []
    wsSF2 = []
    wsNSF = []

    for step in range(0, hklsteps+1):
        hkl = hklbegin + (hklend-hklbegin)/hklsteps*step
        branches = getE(_idxfilehandle, _datafilehandle, hkl)

        newEs = [ branch[0] for branch in branches ]
        Es.extend(newEs)
        wsSF1.extend([ branch[1] if branch[1] > eps else 0. for branch in branches ])
        wsSF2.extend([ branch[2] if branch[2] > eps else 0. for branch in branches ])
        wsNSF.extend([ branch[3] if branch[3] > eps else 0. for branch in branches ])
        qs_h.extend([hkl[0]] * len(newEs))
        qs_k.extend([hkl[1]] * len(newEs))
        qs_l.extend([hkl[2]] * len(newEs))

    # convert to np
    qs_h = np.array(qs_h)
    qs_k = np.array(qs_k)
    qs_l = np.array(qs_l)
    Es = np.array(Es)
    wsSF1 = np.abs(np.array(wsSF1))
    wsSF2 = np.abs(np.array(wsSF2))
    wsNSF = np.abs(np.array(wsNSF))
    wsTotal = wsNSF + wsSF1 + wsSF2

    fig = plot.figure()

    pltSF1 = fig.add_subplot(221)
    pltSF1.set_xlabel("q (rlu)")
    pltSF1.set_ylabel("E (meV)")
    pltSF1.set_xlim(xlim1, xlim2)
    pltSF1.set_ylim(ylim1, ylim2)
    pltSF1.set_title("SF1")
    if plot_dir == 0:
        pltSF1.scatter(qs_h, Es, marker=".", s=wsSF1*symscale)
    elif plot_dir == 1:
        pltSF1.scatter(qs_k, Es, marker=".", s=wsSF1*symscale)
    else:
        pltSF1.scatter(qs_l, Es, marker=".", s=wsSF1*symscale)

    pltSF2 = fig.add_subplot(222)
    pltSF2.set_xlabel("q (rlu)")
    pltSF2.set_ylabel("E (meV)")
    pltSF2.set_xlim(xlim1, xlim2)
    pltSF2.set_ylim(ylim1, ylim2)
    pltSF2.set_title("SF2")
    if plot_dir == 0:
        pltSF2.scatter(qs_h, Es, marker=".", s=wsSF2*symscale)
    elif plot_dir == 1:
        pltSF2.scatter(qs_k, Es, marker=".", s=wsSF2*symscale)
    else:
        pltSF2.scatter(qs_l, Es, marker=".", s=wsSF2*symscale)

    pltNSF = fig.add_subplot(223)
    pltNSF.set_xlabel("q (rlu)")
    pltNSF.set_ylabel("E (meV)")
    pltNSF.set_xlim(xlim1, xlim2)
    pltNSF.set_ylim(ylim1, ylim2)
    pltNSF.set_title("NSF")
    if plot_dir == 0:
        pltNSF.scatter(qs_h, Es, marker=".", s=wsNSF*symscale)
    elif plot_dir == 1:
        pltNSF.scatter(qs_k, Es, marker=".", s=wsNSF*symscale)
    else:
        pltNSF.scatter(qs_l, Es, marker=".", s=wsNSF*symscale)

    pltTotal = fig.add_subplot(224)
    pltTotal.set_xlabel("q (rlu)")
    pltTotal.set_ylabel("E (meV)")
    pltTotal.set_xlim(xlim1, xlim2)
    pltTotal.set_ylim(ylim1, ylim2)
    pltTotal.set_title("Total")
    if plot_dir == 0:
        pltTotal.scatter(qs_h, Es, marker=".", s=wsTotal*symscale)
    elif plot_dir == 1:
        pltTotal.scatter(qs_k, Es, marker=".", s=wsTotal*symscale)
    else:
        pltTotal.scatter(qs_l, Es, marker=".", s=wsTotal*symscale)

    plot.tight_layout()
    plot.show()



# get the energies and weights of the dispersion branches at a specific (h,k,l) coordinate
def get_branch(_idxfilehandle, _datafilehandle, hkl):
    eps = 1e-8

    qs_h = []
    qs_k = []
    qs_l = []
    Es = []
    wsSF1 = []
    wsSF2 = []
    wsNSF = []

    branches = getE(_idxfilehandle, _datafilehandle, hkl)

    newEs = [ branch[0] for branch in branches ]
    Es.extend(newEs)
    wsSF1.extend([ branch[1] if branch[1] > eps else 0. for branch in branches ])
    wsSF2.extend([ branch[2] if branch[2] > eps else 0. for branch in branches ])
    wsNSF.extend([ branch[3] if branch[3] > eps else 0. for branch in branches ])
    qs_h.extend([hkl[0]] * len(newEs))
    qs_k.extend([hkl[1]] * len(newEs))
    qs_l.extend([hkl[2]] * len(newEs))

    # convert to np
    qs_h = np.array(qs_h)
    qs_k = np.array(qs_k)
    qs_l = np.array(qs_l)
    Es = np.array(Es)
    wsSF1 = np.abs(np.array(wsSF1))
    wsSF2 = np.abs(np.array(wsSF2))
    wsNSF = np.abs(np.array(wsNSF))
    wsTotal = wsNSF + wsSF1 + wsSF2

    return [ Es, wsTotal ]



# open index and data file for mapping
try:
    _idxfilehandle = open(idxfile, "rb")
    _datafilehandle = open(datafile, "rb")
except err as IOError:
    print(err)
    exit(-1)



# get the values at a specific coordinate
#print(get_branch(_idxfilehandle, _datafilehandle, [ -0.05, -0.05, 0. ]))
#print(getE(_idxfilehandle, _datafilehandle, [-0.05, -0.05, 0.]))
#exit()



plot_disp(_idxfilehandle, _datafilehandle, plot_hklbegin, plot_hklend, plot_hklsteps)
