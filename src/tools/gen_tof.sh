#!/bin/bash
#
# calculates the polarisation from the tof files
# @author Tobias Weber <tweber@ill.fr>
# @date may-2022
# @license GPLv2 (see 'LICENSE' file)
#

# tools
TOF=./tof

${TOF} \
	../data/mlz_reseda/0p0000941ns/00168499.tof \
	../data/mlz_reseda/0p0000941ns/00168500.tof \
	../data/mlz_reseda/0p0000941ns/00168501.tof \
	../data/mlz_reseda/0p000565ns/00168496.tof \
	../data/mlz_reseda/0p000565ns/00168497.tof \
	../data/mlz_reseda/0p000565ns/00168498.tof \
	../data/mlz_reseda/a_0p00104ns/00168291.tof \
	../data/mlz_reseda/aa_0p00297ns/00168450.tof \
	../data/mlz_reseda/aaa_0p00452ns/00168493.tof \
	../data/mlz_reseda/aaa_0p00452ns/00168494.tof \
	../data/mlz_reseda/aaa_0p00452ns/00168495.tof \
	../data/mlz_reseda/b_0p00971ns/00168292.tof \
	../data/mlz_reseda/bb_0p0122ns/00168434.tof \
	../data/mlz_reseda/bb_0p0122ns/00168435.tof \
	../data/mlz_reseda/bc_0p015ns/00168438.tof \
	../data/mlz_reseda/bc_0p015ns/00168439.tof \
	../data/mlz_reseda/c_0p0204ns/00168360.tof \
	../data/mlz_reseda/c_0p0204ns/00168361.tof \
	../data/mlz_reseda/c_0p0204ns/00168362.tof \
	../data/mlz_reseda/c_0p0204ns/00168363.tof \
	../data/mlz_reseda/c_0p0204ns/00168364.tof \
	../data/mlz_reseda/cc_0p025ns/00168442.tof \
	../data/mlz_reseda/cc_0p025ns/00168449.tof \
	../data/mlz_reseda/ccc_0p0296ns/00168502.tof \
	../data/mlz_reseda/ccc_0p0296ns/00168503.tof \
	../data/mlz_reseda/ccc_0p0296ns/00168504.tof \
	../data/mlz_reseda/d_0p0342ns/00168365.tof \
	../data/mlz_reseda/d_0p0342ns/00168366.tof \
	../data/mlz_reseda/d_0p0342ns/00168367.tof \
	../data/mlz_reseda/d_0p0342ns/00168368.tof \
	../data/mlz_reseda/d_0p0342ns/00168369.tof \
	../data/mlz_reseda/e_0p045ns/00168370.tof \
	../data/mlz_reseda/e_0p045ns/00168371.tof \
	../data/mlz_reseda/e_0p045ns/00168372.tof \
	../data/mlz_reseda/e_0p045ns/00168374.tof \
	../data/mlz_reseda/e_0p045ns/00168375.tof \
	../data/mlz_reseda/ea_0p0512ns/00168505.tof \
	../data/mlz_reseda/ea_0p0512ns/00168507.tof \
	../data/mlz_reseda/ee_0p0574ns/00168490.tof \
	../data/mlz_reseda/ee_0p0574ns/00168491.tof \
	../data/mlz_reseda/ee_0p0574ns/00168492.tof \
	../data/mlz_reseda/f_0p0635ns/00168377.tof \
	../data/mlz_reseda/f_0p0635ns/00168378.tof \
	../data/mlz_reseda/f_0p0635ns/00168379.tof \
	../data/mlz_reseda/f_0p0635ns/00168382.tof \
	../data/mlz_reseda/f_0p0635ns/00168383.tof \
	../data/mlz_reseda/g_0p08ns/00168384.tof \
	../data/mlz_reseda/g_0p08ns/00168385.tof \
	../data/mlz_reseda/g_0p08ns/00168432.tof \
	../data/mlz_reseda/g_0p08ns/00168433.tof \
	../data/mlz_reseda/h_0p1ns/00168293.tof \
	../data/mlz_reseda/i_0p129ns/00168294.tof \
	../data/mlz_reseda/ii_0p146ns/00168508.tof \
	../data/mlz_reseda/ii_0p146ns/00168509.tof \
	../data/mlz_reseda/ii_0p146ns/00168510.tof \
	../data/mlz_reseda/j_0p175ns/00168295.tof \
	../data/mlz_reseda/jj_0p205ns/00168511.tof \
	../data/mlz_reseda/jj_0p205ns/00168512.tof \
	../data/mlz_reseda/k_0p23ns/00168296.tof \
	../data/mlz_reseda/l_0p373ns/00168297.tof \
	../data/mlz_reseda/m_0p542ns/00168298.tof \
	../data/mlz_reseda/m_0p542ns/00168354.tof \
	../data/mlz_reseda/mm_0p694ns/00168443.tof \
	../data/mlz_reseda/mm_0p694ns/00168444.tof \
	../data/mlz_reseda/mm_0p694ns/00168445.tof \
	../data/mlz_reseda/mm_0p694ns/00168446.tof \
	../data/mlz_reseda/mm_0p694ns/00168447.tof \
	../data/mlz_reseda/mm_0p694ns/00168448.tof \
	../data/mlz_reseda/n_0p818ns/00168299.tof \
	../data/mlz_reseda/n_0p818ns/00168356.tof \
	../data/mlz_reseda/n_0p818ns/00168357.tof \
	../data/mlz_reseda/n_0p818ns/00168359.tof \
	../data/mlz_reseda/n_0p818ns/00168488.tof \
	../data/mlz_reseda/nn_1p13ns/00168481.tof \
	../data/mlz_reseda/nn_1p13ns/00168482.tof \
	../data/mlz_reseda/nn_1p13ns/00168483.tof \
	../data/mlz_reseda/nn_1p13ns/00168484.tof \
	../data/mlz_reseda/nn_1p13ns/00168486.tof \
	../data/mlz_reseda/o_1p29ns/00168353.tof \
	../data/mlz_reseda/p_1p61ns/00168300.tof \
	../data/mlz_reseda/q_1p84ns/00168304.tof \
	../data/mlz_reseda/r_2p48ns/00168306.tof \
	../data/mlz_reseda/r_2p48ns/00168307.tof

#	../data/mlz_reseda/aa_0p00297ns/00168451.tof \
#	../data/mlz_reseda/aa_0p00297ns/00168452.tof \
#	../data/mlz_reseda/bb_0p0122ns/00168436.tof \
#	../data/mlz_reseda/bc_0p015ns/00168437.tof \
#	../data/mlz_reseda/cc_0p025ns/00168440.tof \
#	../data/mlz_reseda/cc_0p025ns/00168441.tof \
#	../data/mlz_reseda/e_0p045ns/00168376.tof \
#	../data/mlz_reseda/ea_0p0512ns/00168506.tof \
#	../data/mlz_reseda/jj_0p205ns/00168513.tof \
#	../data/mlz_reseda/n_0p818ns/00168358.tof \
#	../data/mlz_reseda/n_0p818ns/00168487.tof \
#	../data/mlz_reseda/n_0p818ns/00168489.tof \
#	../data/mlz_reseda/nn_1p13ns/00168485.tof \
