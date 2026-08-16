#!/usr/bin/env python3
#
# simple fourier trafo tests
# @author Tobias Weber <tweber@ill.fr>
# @date 16-aug-2026
# @license see 'LICENSE' file
#

import numpy as np
import scipy as sp
import scipy.constants as co
import matplotlib
import matplotlib.pyplot as plt


# settings
eps = 1e-8
sigma_skx = 1e-3
sigma_heli = 1e-3
sigma_elast = 0.5e-3
w_elast = 2e3
min_E = 1e-8
T = 28.5

# domain population
skx_pop = 0.55
heli_pop = 1. - skx_pop


# constants
hbar_in_meVps = co.Planck/co.elementary_charge*1e15/2./np.pi
kB_in_meV_per_K = co.k / co.e * 1e3


# fourier trafo of a gaussian
# see: https://mathworld.wolfram.com/FourierTransformGaussian.html
def gauss_ft(E0, sigma, amp, t):
    return amp * sigma * np.sqrt(2.*np.pi) \
        * np.exp(-0.5 * (t*sigma/hbar_in_meVps)**2.) \
        * np.exp(-1j * t*E0/hbar_in_meVps)


# bose factor
def bose(E, T):
	n = 1. / (np.exp(abs(E)/(kB_in_meV_per_K*T)) - 1.)
	if E >= 0.:
		n += 1.
	return n


dat_skx = np.loadtxt("weightsum_skx.dat")
Es_skx = dat_skx[:, 4]
ws_skx = dat_skx[:, 5] + dat_skx[:, 7]*0.5

dat_heli = np.loadtxt("weightsum_heli.dat")
Es_heli = dat_heli[:, 4]
ws_heli = dat_heli[:, 5] + dat_heli[:, 6] + dat_heli[:, 7]


ts = np.logspace(0.1, 4., 128)
Cs = np.zeros(128)

# elastic peak
Cs += np.real(gauss_ft(0., sigma_elast, w_elast, ts))

# magnons
for E, w in zip(Es_skx, ws_skx):
    if np.abs(E) < min_E or skx_pop < eps:
        continue
    w *= skx_pop
    #w *= bose(E, T)
    Cs += np.real(gauss_ft(E, sigma_skx, w, ts))

for E, w in zip(Es_heli, ws_heli):
    if np.abs(E) < min_E or heli_pop < eps:
        continue
    w *= heli_pop
    #w *= bose(E, T)
    Cs += np.real(gauss_ft(E, sigma_heli, w, ts))

Cs /= Cs[0]


plt.xlabel("t (ps)")
plt.ylabel("C")

plt.semilogx(ts, Cs)
plt.show()
