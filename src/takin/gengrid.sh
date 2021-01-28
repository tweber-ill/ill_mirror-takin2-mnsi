#!/bin/bash
#
# Calculates the dynamical structure factor on a grid
# @author Tobias Weber <tweber@ill.fr>
# @date jan-21
# @license GPLv2 (see 'LICENSE' file)
#

TOOL=./genheli
IDX_START=0
IDX_END=100

for ((idx=$IDX_START; idx<$IDX_END; ++idx)); do
	$TOOL $idx
done
