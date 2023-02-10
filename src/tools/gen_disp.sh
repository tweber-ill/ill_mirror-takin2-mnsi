#!/bin/bash
#
# creates dispersion plots
# @author Tobias Weber <tweber@ill.fr>
# @date 9-feb-2023
# @license GPLv2 (see 'LICENSE' file)
#


# tools
DYN=./dyn
GPL=gnuplot
PLOT_SCRIPT=../src/tools/plot_dyn_path.gpl


dyn_ty=s
T=28.5
B=0.2
num_pts=256


# symmetry points, calculated with symm.jl
Gx=1
Gy=1
Gz=0

K1x=+0.02001
K1y=-0.02001
K1z=+0.0

K2x=+0.01
K2y=-0.01
K2z=-0.0245

Mx=+0.015
My=-0.015
Mz=-0.01225


# G -> M
$DYN --dyntype=$dyn_ty --use_para_perp_calc=0 --outfile=GM.dat \
	--Gx=$Gx --Gy=$Gy --Gz=$Gz \
	--Bx=1 --By=1 --Bz=0 \
	--Px=1 --Py=-1 --Pz=0 \
	--T=$T --B=$B \
	--num_points=$num_pts \
	--qh_start=0 --qk_start=0 --ql_start=0 \
	--qh_end=$Mx --qk_end=$My --ql_end=$Mz


# M -> K
$DYN --dyntype=$dyn_ty --use_para_perp_calc=0 --outfile=MK.dat \
	--Gx=$Gx --Gy=$Gy --Gz=$Gz \
	--Bx=1 --By=1 --Bz=0 \
	--Px=1 --Py=-1 --Pz=0 \
	--T=$T --B=$B \
	--num_points=$num_pts \
	--qh_start=$Mx --qk_start=$My --ql_start=$Mz \
	--qh_end=$K2x --qk_end=$K2y --ql_end=$K2z


# K -> G
$DYN --dyntype=$dyn_ty --use_para_perp_calc=0 --outfile=KG.dat \
	--Gx=$Gx --Gy=$Gy --Gz=$Gz \
	--Bx=1 --By=1 --Bz=0 \
	--Px=1 --Py=-1 --Pz=0 \
	--T=$T --B=$B \
	--num_points=$num_pts \
	--qh_start=$K2x --qk_start=$K2y --ql_start=$K2z \
	--qh_end=0 --qk_end=0 --ql_end=0
