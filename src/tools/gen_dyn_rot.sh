#!/bin/bash
#
# creates dispersion plots rotating from q_para to q_perp
# @author Tobias Weber <tweber@ill.fr>
# @date may-2022
# @license GPLv2 (see 'LICENSE' file)
#

# tools
DYN=./dyn
PY=python3
GPL=gnuplot
PLOT_SCRIPT=../src/tools/plot_dyn_path.gpl
SCANPOS_SCRIPT=../src/tools/draw_recipsetup.gpl
CONV=convert
MPG=ffmpeg

# indices and angles
IDX_START=0     # 0 to 90 degrees
IDX_END=180
IDX_SCALE=0.5   # half-degree steps

calc=1
create_plots=1
create_movie=1


for ((idx=$IDX_START; idx<=$IDX_END; ++idx)); do
	OUTFILE="dyn_${idx}.dat"
	PLOT_FILE="${OUTFILE/.dat/.png}"
	SCANPOS_FILE="${OUTFILE/.dat/_scanpos.png}"
	UNITED_FILE="${OUTFILE/.dat/_united.png}"
	ANGLE=$(${PY} -c "print(${idx} * ${IDX_SCALE})")

	echo -e "\n\x1b[1;34m"
	echo -e "================================================================================"
	echo -e "Calculating ${OUTFILE} with rotation angle ${ANGLE}..."
	echo -e "================================================================================"
	echo -e "\x1b[0m"

	if [ $calc -ne 0 ]; then
		${DYN} --dyntype=h --use_para_perp_calc=0 --outfile="${OUTFILE}" \
			--Gx=1 --Gy=1 --Gz=0 \
			--Bx=1 --By=1 --Bz=0 \
			--Px=1 --Py=-1 --Pz=0 \
			--T=28.5 --B=0.17 \
			--num_points=1024 --explicit_calc=1 \
			--qh_start=-0.1 --qk_start=-0.1 --ql_start=0 \
			--qh_end=0.1 --qk_end=0.1 --ql_end=0 \
			--Rx=0 --Ry=0 --Rz=1 --Ralpha=${ANGLE}
	fi

	if [ $create_plots -ne 0 ]; then
		${GPL} -e "file_dyn = \"${OUTFILE}\"; file_out = \"${PLOT_FILE}\"; out_term = 2;" ${PLOT_SCRIPT}
		${GPL} -e "file_out = \"${SCANPOS_FILE}\"; out_term = 2; scan_angle = ${ANGLE};" ${SCANPOS_SCRIPT}

		echo -e "Uniting ${PLOT_FILE} + ${SCANPOS_FILE} -> ${UNITED_FILE}."
		${CONV} ${PLOT_FILE} "(" ${SCANPOS_FILE} -resize x600 ")" \
			-gravity northeast -geometry +650+120 -composite ${UNITED_FILE}
	fi

	echo -e "\x1b[1;34m"
	echo -e "================================================================================"
	echo -e "\x1b[0m"
done


if [ $create_movie -ne 0 ]; then
	echo -e "\n\x1b[1;34m"
	echo -e "================================================================================"
	echo -e "Creating dispersion video..."
	echo -e "================================================================================"
	echo -e "\x1b[0m"

	${MPG} -i dyn_%d_united.png -y dyn.mp4

	echo -e "\x1b[1;34m"
	echo -e "================================================================================"
	echo -e "\x1b[0m"
fi


# TODO: automatically unite dispersion images
#for ((idx=$IDX_START; idx<=$IDX_END; ++idx)); do
#	let idx_vert=$idx+$IDX_END+1
#
#	printf -v filename_skx "skx/dyn_%d_united.png" ${idx}
#	printf -v filename_heli "heli/dyn_%d_united.png" ${idx}
#	printf -v filename_joined "dyn_joined_%d.png" ${idx}
#
#	printf -v filename_skx_vert "skx_vert/dyn_%d_united.png" ${idx}
#	printf -v filename_heli_vert "heli_vert/dyn_%d_united.png" ${idx}
#	printf -v filename_joined_vert "joined/dyn_joined_%d.png" ${idx_vert}
#
#	if [ ! -f "${filename_skx}" ] || [ ! -f "${filename_heli}" ]; then
#		break
#	fi
#
#	echo -e "Joining ${filename_skx} + ${filename_heli} -> ${filename_joined}"
#	${CONV} +append ${filename_skx} ${filename_heli} ${filename_joined}
#
#	if [ ! -f "${filename_skx_vert}" ] || [ ! -f "${filename_heli_vert}" ]; then
#		break
#	fi
#
#	echo -e "${filename_skx_vert} + ${filename_heli_vert} -> ${filename_joined_vert}"
#	convert +append ${filename_skx_vert} ${filename_heli_vert} ${filename_joined_vert}
#done
