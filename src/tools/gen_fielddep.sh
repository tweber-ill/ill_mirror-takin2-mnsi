#
# calculates the field-dependent magnon energies
# @author Tobias Weber <tweber@ill.fr>
# @date 14-jun-2024
# @license GPLv2 (see 'LICENSE' file)
#

CALC=./weight

Bs_heli=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55)
Bs_fp=(0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5)

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


for B_idx in ${!Bs_heli[@]}; do
	echo "================================================================================"
	echo "B = ${Bs_heli[$B_idx]} T"

		# --outfile="weight_heli_p${B_idx}.dat" \
	${CALC} --dyntype=h --outfile="weight_heli_p${Bs_heli[$B_idx]}.dat" \
		--explicit_calc 1 --do_proj ${proj} \
		--Gx=${Gx} --Gy=${Gy} --Gz=${Gz} \
		--Bx=${Bx} --By=${By} --Bz=${Bz} \
		--Qx=${Qx} --Qy=${Qy} --Qz=${Qz} \
		--T=${T} --B=${Bs_heli[$B_idx]}

		# --outfile="weight_heli_m${B_idx}.dat"
	${CALC} --dyntype=h --outfile="weight_heli_m${Bs_heli[$B_idx]}.dat" \
		--explicit_calc 1 --do_proj ${proj} \
		--Gx=${Gx} --Gy=${Gy} --Gz=${Gz} \
		--Bx=-${Bx} --By=-${By} --Bz=-${Bz} \
		--Qx=${Qx} --Qy=${Qy} --Qz=${Qz} \
		--T=${T} --B=${Bs_heli[$B_idx]}
	echo "================================================================================\n\n"
done


for B_idx in ${!Bs_fp[@]}; do
	echo "================================================================================"
	echo "B = ${Bs_fp[$B_idx]} T"

		# --outfile="weight_fp_p${B_idx}.dat"
	${CALC} --dyntype=f --outfile="weight_fp_p${Bs_fp[$B_idx]}.dat" \
		--explicit_calc 1 --do_proj ${proj} \
		--Gx=${Gx} --Gy=${Gy} --Gz=${Gz} \
		--Bx=${Bx} --By=${By} --Bz=${Bz} \
		--Qx=${Qx} --Qy=${Qy} --Qz=${Qz} \
		--T=${T} --B=${Bs_fp[$B_idx]}

		# --outfile="weight_fp_m${B_idx}.dat"
	${CALC} --dyntype=f --outfile="weight_fp_m${Bs_fp[$B_idx]}.dat" \
		--explicit_calc 1 --do_proj ${proj} \
		--Gx=${Gx} --Gy=${Gy} --Gz=${Gz} \
		--Bx=-${Bx} --By=-${By} --Bz=-${Bz} \
		--Qx=${Qx} --Qy=${Qy} --Qz=${Qz} \
		--T=${T} --B=${Bs_fp[$B_idx]}
	echo "================================================================================\n\n"
done


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
echo "================================================================================\n\n"
