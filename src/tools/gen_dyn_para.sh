#!/bin/bash
#
# creates dispersion plots comparing skyrmion modes with conical band formation when increasing q_perp
# @author Tobias Weber <tweber@ill.fr>
# @date may-2022
# @license GPLv2 (see 'LICENSE' file)
#

# tools
DYN=./dyn
PY=python3
GPL=gnuplot
PLOT_SCRIPT=../src/tools/plot_dyn.gpl
MPG=ffmpeg

# indices and angles
IDX_START=0
IDX_END=150
IDX_SCALE=0.0005

calc=1
create_plots=1
create_movie=1

along_qpara=1   # flip q_para und q_perp directions


for ((idx=$IDX_START; idx<=$IDX_END; ++idx)); do
	QPERP=$(${PY} -c "print(\"%.4f\" % (${idx} * ${IDX_SCALE}))")

	OUTFILE="dyn_${idx}.dat"
	PLOTFILE="${OUTFILE/.dat/.png}"

	#OUTFILE_SKX="dyn_skx_${idx}.dat"
	#PLOTFILE_SKX="${OUTFILE_SKX/.dat/.png}"

	echo -e "\n\x1b[1;34m"
	echo -e "================================================================================"
	echo -e "Calculating ${OUTFILE} with q_perp = ${QPERP}..."
	echo -e "================================================================================"
	echo -e "\x1b[0m"

	if [ $calc -ne 0 ]; then
		${DYN} --dyntype=h --use_para_perp_calc=1 --outfile="${OUTFILE}" \
			--Gx=1 --Gy=1 --Gz=0 \
			--Bx=1 --By=1 --Bz=0 \
			--Px=1 --Py=-1 --Pz=0 \
			--T=28.5 --B=0.15 \
			--explicit_calc=1 --along_qpara=${along_qpara} \
			--qperpx=1 --qperpy=-1 --qperpz=0 --qperp=${QPERP} \
			--qrange=0.12 --qdelta=0.00024

		#${DYN} --dyntype=s --use_para_perp_calc=1 --outfile="${OUTFILE_SKX}" \
		#	--Gx=1 --Gy=1 --Gz=0 \
		#	--Bx=1 --By=1 --Bz=0 \
		#	--Px=1 --Py=-1 --Pz=0 \
		#	--T=28.5 --B=0.15 \
		#	--explicit_calc=1 --along_qpara=${along_qpara} \
		#	--qperpx=1 --qperpy=-1 --qperpz=0 --qperp=${QPERP} \
		#	--qrange=0.12 --qdelta=0.00024
	fi

	if [ $create_plots -ne 0 ]; then
		${GPL} -e "file_dyn = \"${OUTFILE}\"; file_out = \"${PLOTFILE}\"; out_term = 2; along_q_para=\"${along_qpara}\"; q_perp = \"${QPERP}\";" ${PLOT_SCRIPT}
		#${GPL} -e "file_dyn = \"${OUTFILE_SKX}\"; file_out = \"${PLOTFILE_SKX}\"; out_term = 2; along_q_para=\"${along_qpara}\"; q_perp = \"${QPERP}\";" ${PLOT_SCRIPT}
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

	${MPG} -framerate 20 -i dyn_%d.png -y dyn.mp4
	#${MPG} -framerate 20 -i dyn_skx_%d.png -y dyn_skx.mp4

	echo -e "\x1b[1;34m"
	echo -e "================================================================================"
	echo -e "\x1b[0m"
fi
