#!/bin/bash
#
# creates dispersion plots rotating around an axis
# @author Tobias Weber <tweber@ill.fr>
# @date may-2022
# @license GPLv2 (see 'LICENSE' file)
#

DYN=./dyn
PY=python

IDX_START=0    # 0 to 90 degrees
IDX_END=180
IDX_SCALE=0.5  # half-angle steps

for ((idx=$IDX_START; idx<$IDX_END; ++idx)); do
	OUTFILE="dyn_${idx}.dat"
	ANGLE=$(${PY} -c "print(${idx} * ${IDX_SCALE})")

	echo -e "\n"
	echo -e "--------------------------------------------------------------------------------"
	echo -e "Calculating ${OUTFILE} with rotation angle ${ANGLE}..."
	echo -e "--------------------------------------------------------------------------------"

	${DYN} --dyntype=h --use_para_perp_calc=0 --outfile="${OUTFILE}" \
		--Gx=1 --Gy=1 --Gz=0 \
		--Bx=1 --By=1 --Bz=0 \
		--Px=1 --Py=-1 --Pz=0 \
		--T=28.5 --B=0.15 \
		--num_points=256 \
		--qh_start=-0.2 --qk_start=-0.2 --ql_start=0 \
		--qh_end=0.2 --qk_end=0.2 --ql_end=0 \
		--Rx=0 --Ry=0 --Rz=1 --Ralpha=${ANGLE}

	echo -e "--------------------------------------------------------------------------------\n"
done
