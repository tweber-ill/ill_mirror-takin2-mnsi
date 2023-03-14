#
# calculates the multi-domain helical dispersion
# @author Tobias Weber <tweber@ill.fr>
# @date mar-2023
# @license GPLv2 (see 'LICENSE' file)
#

calc_disp=1
plot_disp=1


DYN=./dyn
GPL=gnuplot
PLOT_SCRIPT=../src/tools/plot_multidomain.sh

qh_start=-0.06
qk_start=-0.06
ql_start=0

qh_end=0.06
qk_end=0.06
ql_end=0

T=20
B_mag=0.001


if [ $calc_disp -ne 0 ]; then
	${DYN} --dyntype=h --explicit_calc 1 \
		--use_para_perp_calc=0  --do_proj 1 \
		--outfile=heli_domain1.dat \
		--Gx=1 --Gy=1 --Gz=0 \
		--Bx=1 --By=1 --Bz=1 --T=${T} --B=${B_mag} \
		--num_points=512 \
		--qh_start=${qh_start} --qk_start=${qk_start} --ql_start=${ql_start} \
		--qh_end=${qh_end} --qk_end=${qk_end} --ql_end=${ql_end}

	${DYN} --dyntype=h --explicit_calc 1 \
		--use_para_perp_calc=0  --do_proj 1 \
		--outfile=heli_domain2.dat \
		--Gx=1 --Gy=1 --Gz=0 \
		--Bx=-1 --By=1 --Bz=1 --T=${T} --B=${B_mag} \
		--num_points=512 \
		--qh_start=${qh_start} --qk_start=${qk_start} --ql_start=${ql_start} \
		--qh_end=${qh_end} --qk_end=${qk_end} --ql_end=${ql_end}

	${DYN} --dyntype=h --explicit_calc 1 \
		--use_para_perp_calc=0  --do_proj 1 \
		--outfile=heli_domain3.dat \
		--Gx=1 --Gy=1 --Gz=0 \
		--Bx=1 --By=-1 --Bz=1 --T=${T} --B=${B_mag} \
		--num_points=512 \
		--qh_start=${qh_start} --qk_start=${qk_start} --ql_start=${ql_start} \
		--qh_end=${qh_end} --qk_end=${qk_end} --ql_end=${ql_end}

	${DYN} --dyntype=h --explicit_calc 1 \
		--use_para_perp_calc=0  --do_proj 1 \
		--outfile=heli_domain4.dat \
		--Gx=1 --Gy=1 --Gz=0 \
		--Bx=1 --By=1 --Bz=-1 --T=${T} --B=${B_mag} \
		--num_points=512 \
		--qh_start=${qh_start} --qk_start=${qk_start} --ql_start=${ql_start} \
		--qh_end=${qh_end} --qk_end=${qk_end} --ql_end=${ql_end}
fi


if [ $plot_disp -ne 0 ]; then
	${GPL} -e \
		"file_domain1 = \"heli_domain1.dat\"; \
		 file_domain2 = \"heli_domain2.dat\"; \
		 file_domain3 = \"heli_domain3.dat\"; \
		 file_domain4 = \"heli_domain4.dat\"; \
		file_out = \"heli_multidomain.pdf\";" \
		${PLOT_SCRIPT}
fi
