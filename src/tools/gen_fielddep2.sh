#!/bin/bash
#
# calculates the field-dependent magnon energies
# @author Tobias Weber <tweber@ill.fr>
# @date 14-jun-2024
# @license GPLv2 (see 'LICENSE' file)
#

CALC=./weight
CALC_B=./fielddep

Gx=1
Gy=1
Gz=0

Bx=1
By=0
Bz=0

Qx=1.0
Qy=1.0001
Qz=0.0

proj=0
T=5

B_start_fp=0.58
B_end_fp=1.5
B_start_heli=0
B_end_heli=0.58
B_step=0.001


echo "================================================================================"
${CALC_B} --dyntype=f --outfile="weight_fp_p.dat" \
	--explicit_calc 1 --do_proj ${proj} \
	--Gx=${Gx} --Gy=${Gy} --Gz=${Gz} \
	--Bx=${Bx} --By=${By} --Bz=${Bz} \
	--Qx=${Qx} --Qy=${Qy} --Qz=${Qz} \
	--T=${T} \
	--B_start=${B_start_fp} --B_end=${B_end_fp} --B_step=${B_step}

${CALC_B} --dyntype=f --outfile="weight_fp_m.dat" \
	--explicit_calc 1 --do_proj ${proj} \
	--Gx=${Gx} --Gy=${Gy} --Gz=${Gz} \
	--Bx=-${Bx} --By=-${By} --Bz=-${Bz} \
	--Qx=${Qx} --Qy=${Qy} --Qz=${Qz} \
	--T=${T} \
	--B_start=${B_start_fp} --B_end=${B_end_fp} --B_step=${B_step}
echo "================================================================================"



echo "================================================================================"
${CALC_B} --dyntype=h --outfile="weight_heli_p.dat" \
	--explicit_calc 1 --do_proj ${proj} \
	--Gx=${Gx} --Gy=${Gy} --Gz=${Gz} \
	--Bx=${Bx} --By=${By} --Bz=${Bz} \
	--Qx=${Qx} --Qy=${Qy} --Qz=${Qz} \
	--T=${T} \
	--B_start=${B_start_heli} --B_end=${B_end_heli} --B_step=${B_step}

${CALC_B} --dyntype=h --outfile="weight_heli_m.dat" \
	--explicit_calc 1 --do_proj ${proj} \
	--Gx=${Gx} --Gy=${Gy} --Gz=${Gz} \
	--Bx=-${Bx} --By=-${By} --Bz=-${Bz} \
	--Qx=${Qx} --Qy=${Qy} --Qz=${Qz} \
	--T=${T} \
	--B_start=${B_start_heli} --B_end=${B_end_heli} --B_step=${B_step}
echo "================================================================================"



echo "================================================================================"
${CALC} --dyntype=s --outfile="weight_skx_m.dat" \
	--explicit_calc 1 --do_proj ${proj} \
	--Gx=${Gx} --Gy=${Gy} --Gz=${Gz} \
	--Bx=-${Bx} --By=-${By} --Bz=-${Bz} \
	--Qx=${Qx} --Qy=${Qy} --Qz=${Qz} \
	--Px=1 --Py=-1 --Pz=0 \
	--T=${T} --B=0.2

${CALC} --dyntype=s --outfile="weight_skx_p.dat" \
	--explicit_calc 1 --do_proj ${proj} \
	--Gx=${Gx} --Gy=${Gy} --Gz=${Gz} \
	--Bx=${Bx} --By=${By} --Bz=${Bz} \
	--Qx=${Qx} --Qy=${Qy} --Qz=${Qz} \
	--Px=1 --Py=-1 --Pz=0 \
	--T=${T} --B=0.2
echo "================================================================================"
