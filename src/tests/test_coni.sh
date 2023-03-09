#
# compares explicit and implicit conical calculations
# @author Tobias Weber <tweber@ill.fr>
# @date mar-2023
# @license GPLv2 (see 'LICENSE' file)
#

calc_expl=1
calc_impl=1
calc_gs=1


if [ $calc_expl -ne 0 ]; then
	./dyn --dyntype=h --explicit_calc 1 \
		--use_para_perp_calc=1  --along_qpara 1 --do_proj 1 \
		--outfile=dyn_para_expl.dat \
		--Gx=1 --Gy=1 --Gz=0 \
		--Bx=1 --By=1 --Bz=0 --T=20 --B=0.15 \
		--qrange 0.15 --qdelta 0.0005 \
		--qperpx 1 --qperpy -1 --qperpz 0 \
		--qperp 0 --qperp2 0

	gnuplot -e "file_dyn = \"dyn_para_expl.dat\"; file_out = \"dyn_para_expl.pdf\";" \
		../src/tools/plot_dyn.gpl


	./dyn --dyntype=h --explicit_calc 1 \
		--use_para_perp_calc=1  --along_qpara 0 --do_proj 1 \
		--outfile=dyn_perp_expl.dat \
		--Gx=1 --Gy=1 --Gz=0 \
		--Bx=1 --By=1 --Bz=0 --T=20 --B=0.15 \
		--qrange 0.15 --qdelta 0.0005 \
		--qperpx 1 --qperpy -1 --qperpz 0 \
		--qperp 0 --qperp2 0

	gnuplot -e "file_dyn = \"dyn_perp_expl.dat\"; file_out = \"dyn_perp_expl.pdf\";" \
		../src/tools/plot_dyn.gpl
fi


if [ $calc_gs -ne 0 ]; then
	# move to lower B step-by-step to keep the F_min fit valid
	./heli_gs --T_exp 20 --T_theo -1000 --B_theo 22.5 --outfile coni_gs.bin
	./heli_gs --T_exp 20 --T_theo -1000 --B_theo 20 --gsfile coni_gs.bin --outfile coni_gs.bin
	./heli_gs --T_exp 20 --T_theo -1000 --B_theo 18 --gsfile coni_gs.bin --outfile coni_gs.bin
	./heli_gs --T_exp 20 --T_theo -1000 --B_theo 16 --gsfile coni_gs.bin --outfile coni_gs.bin
	./heli_gs --T_exp 20 --T_theo -1000 --B_theo 14 --gsfile coni_gs.bin --outfile coni_gs.bin
fi


if [ $calc_impl -ne 0 ]; then
	./dyn --dyntype=h --explicit_calc 0 --gsfile=coni_gs.bin \
		--use_para_perp_calc=1  --along_qpara 1 --do_proj 1 \
		--outfile=dyn_para_impl.dat \
		--Gx=1 --Gy=1 --Gz=0 \
		--Bx=1 --By=1 --Bz=0 --T=20 --B=0.15 \
		--qrange 0.15 --qdelta 0.0005 \
		--qperpx 1 --qperpy -1 --qperpz 0 \
		--qperp 0 --qperp2 0

	gnuplot -e "file_dyn = \"dyn_para_impl.dat\"; file_out = \"dyn_para_impl.pdf\";" \
		../src/tools/plot_dyn.gpl


	./dyn --dyntype=h --explicit_calc 0 --gsfile=coni_gs.bin \
		--use_para_perp_calc=1  --along_qpara 0 --do_proj 1 \
		--outfile=dyn_perp_impl.dat \
		--Gx=1 --Gy=1 --Gz=0 \
		--Bx=1 --By=1 --Bz=0 --T=20 --B=0.15 \
		--qrange 0.15 --qdelta 0.0005 \
		--qperpx 1 --qperpy -1 --qperpz 0 \
		--qperp 0 --qperp2 0

	gnuplot -e "file_dyn = \"dyn_perp_impl.dat\"; file_out = \"dyn_perp_impl.pdf\";" \
		../src/tools/plot_dyn.gpl
fi
