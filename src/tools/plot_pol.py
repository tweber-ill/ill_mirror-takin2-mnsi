#
# plots polarisations vs. spin-echo time
# @author Tobias Weber <tweber@ill.fr>
# @date may-2022
# @license GPLv2 (see 'LICENSE' file)
#

import numpy as np
import matplotlib.pyplot as plt

P1 = np.genfromtxt("P_1.dat")
P2 = np.genfromtxt("P_2.dat")
P3 = np.genfromtxt("P_3.dat")
P4 = np.genfromtxt("P_4.dat")
P5 = np.genfromtxt("P_5.dat")
P6 = np.genfromtxt("P_6.dat")

P0_1 = np.genfromtxt("P0_1.dat")
P0_2 = np.genfromtxt("P0_2.dat")
P0_3 = np.genfromtxt("P0_3.dat")
P0_4 = np.genfromtxt("P0_4.dat")
P0_5 = np.genfromtxt("P0_5.dat")
P0_6 = np.genfromtxt("P0_6.dat")

taus = P1[:,3]
pols = (
#	  P1[:,4] / P0_1[:,4]
#	+ P2[:,4] / P0_2[:,4]
	+ P3[:,4] / P0_3[:,4]
	+ P4[:,4] / P0_4[:,4]
#	+ P5[:,4] / P0_5[:,4]
#	+ P6[:,4] / P0_6[:,4]
	) / 2.
errs = (
#	  np.sqrt((P1[:,5]/P0_1[:,4])**2. + (P1[:,4]*P0_1[:,5]/P0_1[:,4]**2.)**2.)
#	+ np.sqrt((P2[:,5]/P0_2[:,4])**2. + (P2[:,4]*P0_2[:,5]/P0_2[:,4]**2.)**2.)
	+ np.sqrt((P3[:,5]/P0_3[:,4])**2. + (P3[:,4]*P0_3[:,5]/P0_3[:,4]**2.)**2.)
	+ np.sqrt((P4[:,5]/P0_4[:,4])**2. + (P4[:,4]*P0_4[:,5]/P0_4[:,4]**2.)**2.)
#	+ np.sqrt((P5[:,5]/P0_5[:,4])**2. + (P5[:,4]*P0_5[:,5]/P0_5[:,4]**2.)**2.)
#	+ np.sqrt((P6[:,5]/P0_6[:,4])**2. + (P6[:,4]*P0_6[:,5]/P0_6[:,4]**2.)**2.)
	) / 2.

plt.figure()
plt.xlabel("Spin-echo time τ (ps)")
plt.ylabel("Polarisation")
plt.xscale("log")
plt.ylim([0., 1.])
plt.errorbar(taus, pols, marker="o", yerr=errs)
plt.show()
