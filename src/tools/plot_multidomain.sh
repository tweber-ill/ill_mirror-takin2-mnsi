#!/usr/local/bin/gnuplot -p

#
# plots a multi-domain dispersion
# @author Tobias Weber <tweber@ill.fr>
# @date august-2018
# @license GPLv2 (see 'LICENSE' file)
#

qrange        = 0.06  # q range to plot
Erange        = 0.2   # E range to plot
scale_weight  = 1.5   # weight scale for plot
key_pointsize = 2.0   # symbol size in legend

show_pol      = 0
show_key      = 0     # show legend
show_scans    = 1     # indicate scans


a      = 4.558        # lattice constant
kh_A   = 0.036        # at 20 K
#kh_A   = 0.039        # at 29 K
kh_rlu = kh_A / (2*pi / a)

Bc2    = 0.53
#Bc2    = 0.3223
g      = 2.
muB    = 5.78838e-02  # codata value


if(!exists("out_term")) {
	out_term = 1          # 1: pdf, 2: png, others: screen
}

if(out_term == 1) {
	set term pdf color enhanced font "NimbusSans-Regular, 50" size 20, 15
}
if(out_term == 2) {
	set term pngcairo enhanced font "NimbusSans-Regular, 56" size 2666, 2000

	scale_weight = 2.0 * scale_weight
	key_pointsize = 2.0 * key_pointsize
}

set clip points
set border lw 8


# input data files
if(!exists("file_domain1")) {
	file_domain1 = "domain1.dat"
}
if(!exists("file_domain2")) {
	file_domain1 = "domain2.dat"
}
if(!exists("file_domain3")) {
	file_domain1 = "domain3.dat"
}
if(!exists("file_domain4")) {
	file_domain1 = "domain4.dat"
}

# output plot file
if(!exists("file_out")) {
	if(out_term == 1) {
		file_out = "dyn.pdf"
	}
	if(out_term == 2) {
		file_out = "dyn.png"
	}
}

if(exists("file_out")) {
	set output file_out
	print(sprintf("Plotting %s.", file_out))
}


# colours
col_axis    = "#000000"  # colour for coordinate cross
col_label   = "#000000"  # colour for label
col_scan    = "#aaaaaa"  # colour for scan marker
col_full    = "#999999"  # colour for all weights
col_pol1    = "#0000bb"  # colour for spin-flip channel 1
col_pol2    = "#bb0000"  # colour for spin-flip channel 2
col_pol3    = "#009900"  # colour for non-spin-flip channel
col_domain1 = "#000000"  # colour for domain 1
col_domain2 = "#ff0000"  # colour for domain 1
col_domain3 = "#00aa00"  # colour for domain 1
col_domain4 = "#0000ff"  # colour for domain 1


# ranges
set xrange  [ -qrange               : qrange               ]
set x2range [ -qrange/kh_rlu        : qrange/kh_rlu        ]
set yrange  [ -Erange               : Erange               ]
set y2range [ -Erange / (g*muB*Bc2) : Erange / (g*muB*Bc2) ]


# ticks
set xtics 0.05 nomirror
set x2tics 1
set mxtics 5
set mx2tics 2

set ytics 0.1 nomirror
set y2tics 2
set mytics 5
set my2tics 4


# labels
set xlabel "{/NimbusSans-Italic q} (rlu)" \
	offset 0, 0.25
set x2label "{/NimbusSans-Italic q} (k_h)" \
	offset 0, -0.4

set ylabel \
	"{/NimbusSans-Italic E} (meV)" \
	offset 2, 0
set y2label \
	"{/NimbusSans-Italic E} ({/NimbusSans-Italic g} {/NimbusSans-Italic μ_0} {/NimbusSans-Italic μ_B} {/NimbusSans-Italic H}_{c2}^{int})" \
	offset -0.25, 0


# legend
if(show_key != 0) {
        set key at graph 0.49, graph 0.97 height 0.3 samplen 0 Left reverse box opaque
}
else {
        unset key
}


# show coordinate cross
set arrow 1 from 0,-Erange to 0,Erange lw 2 dt 2 lc rgb col_axis nohead
set arrow 2 from -qrange,0 to qrange,0 lw 2 dt 2 lc rgb col_axis nohead


# indicate scans
if(show_scans != 0) {
	scanw = 0.0025

	qpos1 = -0.04
	qpos2 = +0.04

	E_begin1 = -0.15
	E_end1   = +0.15
	E_begin2 = -0.15
	E_end2   = +0.15

	set obj 1 rect \
		from first qpos1-scanw, first E_begin1 \
		to first qpos1+scanw, first E_end1 \
		fc rgb col_scan fs transparent solid 0.5 front

	set obj 2 rect \
		from first qpos2-scanw, first E_begin2 \
		to first qpos2+scanw, first E_end2 \
		fc rgb col_scan fs transparent solid 0.5 front

	set label 1 at qpos1, 0.075 center "(i)" tc rgb col_label front
	set label 2 at qpos2, 0.075 center "(ii)" tc rgb col_label front
}


# weight factor
weight_raw(w)         = (w > 5.) ? 5. : sqrt(w) * 2.5
weight(w)             = (weight_raw(w) < 0.04) ? 0 : weight_raw(w)
weight_ifminusE(w, E) = (E < 0) ? weight(w) : 0


E_sign = 1.

$dummy_data << END_DATA
	NaN NaN
END_DATA

if(show_pol != 0) {
	# polarised results
	plot \
		file_domain1 using ($9):(E_sign*$4) \
			with points pt 0        lc rgb col_full notitle, \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol2 title "SF (+-)", \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol1 title "SF (-+)", \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol3 title "NSF", \
		file_domain1 using ($9):(E_sign*$4):(weight($7)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol2 notitle, \
		file_domain1 using ($9):(E_sign*$4):(weight($6)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol1 notitle, \
		file_domain1 using ($9):(E_sign*$4):(weight_ifminusE($7, $4)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol2 notitle, \
		file_domain1 using ($9):(E_sign*$4):(weight($8)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol3 notitle, \
		file_domain2 using ($9):(E_sign*$4) \
			with points pt 0        lc rgb col_full notitle, \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol2 notitle, \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol1 notitle, \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol3 notitle, \
		file_domain2 using ($9):(E_sign*$4):(weight($7)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol2 notitle, \
		file_domain2 using ($9):(E_sign*$4):(weight($6)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol1 notitle, \
		file_domain2 using ($9):(E_sign*$4):(weight_ifminusE($7, $4)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol2 notitle, \
		file_domain2 using ($9):(E_sign*$4):(weight($8)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol3 notitle, \
		file_domain3 using ($9):(E_sign*$4) \
			with points pt 0        lc rgb col_full notitle, \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol2 notitle, \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol1 notitle, \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol3 notitle, \
		file_domain3 using ($9):(E_sign*$4):(weight($7)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol2 notitle, \
		file_domain3 using ($9):(E_sign*$4):(weight($6)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol1 notitle, \
		file_domain3 using ($9):(E_sign*$4):(weight_ifminusE($7, $4)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol2 notitle, \
		file_domain3 using ($9):(E_sign*$4):(weight($8)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol3 notitle, \
		file_domain4 using ($9):(E_sign*$4) \
			with points pt 0        lc rgb col_full notitle, \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol2 notitle, \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol1 notitle, \
		"$dummy_data" with points pt 7 ps key_pointsize lc rgb col_pol3 notitle, \
		file_domain4 using ($9):(E_sign*$4):(weight($7)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol2 notitle, \
		file_domain4 using ($9):(E_sign*$4):(weight($6)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol1 notitle, \
		file_domain4 using ($9):(E_sign*$4):(weight_ifminusE($7, $4)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol2 notitle, \
		file_domain4 using ($9):(E_sign*$4):(weight($8)*scale_weight) \
			with points pt 7 ps var lc rgb col_pol3 notitle
}
else {
	# unpolarised results
	plot \
	        file_domain1 using ($9):(E_sign*$4) \
			with points pt 0        lc rgb col_full notitle, \
	        "$dummy_data" with points pt 7 ps key_pointsize lc rgb col_domain1 title "Domain 1", \
	        file_domain1 using ($9):(E_sign*$4):(weight($6+$7+$8)*scale_weight) \
			with points pt 7 ps var lc rgb col_domain1 notitle, \
	        file_domain2 using ($9):(E_sign*$4) \
			with points pt 0        lc rgb col_full notitle, \
	        "$dummy_data" with points pt 7 ps key_pointsize lc rgb col_domain2 title "Domain 2", \
	        file_domain2 using ($9):(E_sign*$4):(weight($6+$7+$8)*scale_weight) \
			with points pt 7 ps var lc rgb col_domain2 notitle, \
	        file_domain3 using ($9):(E_sign*$4) \
			with points pt 0        lc rgb col_full notitle, \
	        "$dummy_data" with points pt 7 ps key_pointsize lc rgb col_domain3 title "Domain 3", \
	        file_domain3 using ($9):(E_sign*$4):(weight($6+$7+$8)*scale_weight) \
			with points pt 7 ps var lc rgb col_domain3 notitle, \
	        file_domain4 using ($9):(E_sign*$4) \
			with points pt 0        lc rgb col_full notitle, \
	        "$dummy_data" with points pt 7 ps key_pointsize lc rgb col_domain4 title "Domain 4", \
	        file_domain4 using ($9):(E_sign*$4):(weight($6+$7+$8)*scale_weight) \
			with points pt 7 ps var lc rgb col_domain4 notitle
}
