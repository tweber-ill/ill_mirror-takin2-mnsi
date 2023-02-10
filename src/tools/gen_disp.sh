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
Gh=1
Gk=1
Gl=0

Mh=+0.01
Mk=-0.01
Ml=+0.0

Kh=+0.01
Kk=-0.01
Kl=-0.00817


# G -> M
$DYN --dyntype=$dyn_ty --use_para_perp_calc=0 --outfile=GM.dat \
	--Gx=$Gh --Gy=$Gk --Gz=$Gl \
	--Bx=1 --By=1 --Bz=0 \
	--Px=1 --Py=-1 --Pz=0 \
	--T=$T --B=$B \
	--num_points=$num_pts \
	--qh_start=0 --qk_start=0 --ql_start=0 \
	--qh_end=$Mh --qk_end=$Mk --ql_end=$Ml


# M -> K
$DYN --dyntype=$dyn_ty --use_para_perp_calc=0 --outfile=MK.dat \
	--Gx=$Gh --Gy=$Gk --Gz=$Gl \
	--Bx=1 --By=1 --Bz=0 \
	--Px=1 --Py=-1 --Pz=0 \
	--T=$T --B=$B \
	--num_points=$num_pts \
	--qh_start=$Mh --qk_start=$Mk --ql_start=$Ml \
	--qh_end=$Kh --qk_end=$Kk --ql_end=$Kl


# K -> G
$DYN --dyntype=$dyn_ty --use_para_perp_calc=0 --outfile=KG.dat \
	--Gx=$Gh --Gy=$Gk --Gz=$Gl \
	--Bx=1 --By=1 --Bz=0 \
	--Px=1 --Py=-1 --Pz=0 \
	--T=$T --B=$B \
	--num_points=$num_pts \
	--qh_start=$Kh --qk_start=$Kk --ql_start=$Kl \
	--qh_end=0 --qk_end=0 --ql_end=0
