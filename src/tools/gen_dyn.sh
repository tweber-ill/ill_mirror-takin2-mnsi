#!/bin/bash
#
# creates dispersion plots rotating around an axis
# @author Tobias Weber <tweber@ill.fr>
# @date may-2022
# @license GPLv2 (see 'LICENSE' file)
#

# tools
DYN=./dyn
PY=python
GPL=gnuplot
PLOT_SCRIPT=../src/tools/plot_dyn_path.gpl
MPG=ffmpeg

# indices and angles
IDX_START=0     # 0 to 90 degrees
IDX_END=180
IDX_SCALE=0.5   # half-degree steps

create_plots=1
create_movie=1


for ((idx=$IDX_START; idx<=$IDX_END; ++idx)); do
	OUTFILE="dyn_${idx}.dat"
	PLOTFILE="${OUTFILE/.dat/.png}"
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
		--num_points=512 \
		--qh_start=-0.2 --qk_start=-0.2 --ql_start=0 \
		--qh_end=0.2 --qk_end=0.2 --ql_end=0 \
		--Rx=0 --Ry=0 --Rz=1 --Ralpha=${ANGLE}

	if [ $create_plots -ne 0 ]; then
		${GPL} -e "file_dyn = \"${OUTFILE}\"; file_out = \"${PLOTFILE}\"; out_term = 2;" ${PLOT_SCRIPT}
	fi

	echo -e "--------------------------------------------------------------------------------\n"
done


if [ $create_movie -ne 0 ]; then
	echo -e "\n"
	echo -e "--------------------------------------------------------------------------------"
	echo -e "Creating dispersion movie..."
	echo -e "--------------------------------------------------------------------------------"

	${MPG} -i dyn_%d.png -y dyn.mp4

	echo -e "--------------------------------------------------------------------------------"
fi
