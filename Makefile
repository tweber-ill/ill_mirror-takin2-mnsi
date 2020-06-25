#
# MnSi dynamics module for Takin
# @author Tobias Weber <tweber@ill.fr>
# @date 2018-2020
# @license GPLv2 (see 'LICENSE' file)
#

CXX = g++
STD = -std=c++17
OPT = -O2 -march=native
DEFS = -DDEF_SKX_ORDER=7 -DDEF_HELI_ORDER=7 \
	-DNO_MINIMISATION -DNO_REDEFINITIONS \
	-D__HACK_FULL_INST__
INCS = -Isrc -Iext -Iext/takin -I/usr/local/include \
	-I/usr/include/lapacke -I/usr/local/opt/lapack/include \
	-I/usr/include/qt5 -I/usr/include/x86_64-linux-gnu/qt5/ \
	#-I/home/tw/build/boost_1_73_0
LIBDIRS = -L/usr/local/opt/lapack/lib -L/usr/local/lib
LIBDEFS = -fPIC


# -----------------------------------------------------------------------------
.PHONY: all clean

all: prepare lib/skxmod.so lib/skxmod_grid.so \
	bin/genskx bin/merge bin/convert bin/dump \
	bin/heliphase bin/skx_gs \
	bin/drawskx bin/dyn bin/weight #bin/weight_sum

clean:
	find . -name "*.o" -exec rm -fv {} \;
	rm -rfv bin/
	rm -rfv lib/

prepare:
	mkdir -p bin/
	mkdir -p lib/
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Takin plugin modules
# -----------------------------------------------------------------------------
lib/skxmod.so: src/core/skx.o src/core/fp.o src/core/heli.o src/core/magsys.o src/takin/takin.o \
	ext/takin/tools/monteconvo/sqwbase.o ext/tlibs2/libs/log.o
	@echo "Linking Takin module $@..."
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -shared -o $@ $+ -llapacke
	strip $@

lib/skxmod_grid.so: src/takin/takin_grid.o ext/takin/tools/monteconvo/sqwbase.o ext/tlibs2/libs/log.o
	@echo "Linking Takin grid module $@..."
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -shared -o $@ $+ -lboost_system -lQt5Core
	strip $@
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# tools
# -----------------------------------------------------------------------------
bin/genskx: src/takin/genskx.o src/core/skx.o src/core/magsys.o ext/tlibs2/libs/log.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -o $@ $+ -lboost_filesystem -llapacke -lpthread
	strip $@

bin/merge: src/takin/merge.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+ -lboost_system
	strip $@

bin/convert: src/takin/convert.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+
	strip $@

bin/dump: src/takin/dump.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+
	strip $@

bin/drawskx: src/calc/drawskx.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+
	strip $@

bin/dyn: src/calc/dyn.o src/core/skx.o src/core/fp.o src/core/heli.o src/core/magsys.o ext/tlibs2/libs/log.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -o $@ $+ -llapacke -lpthread
	strip $@

bin/weight: src/calc/weight.o src/core/skx.o src/core/fp.o src/core/heli.o src/core/magsys.o ext/tlibs2/libs/log.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -o $@ $+ -llapacke -lpthread
	strip $@

bin/weight_sum: src/calc/weight_sum.o src/core/skx.o src/core/fp.o src/core/heli.o src/core/magsys.o ext/tlibs2/libs/log.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -o $@ $+ -llapacke -lpthread
	strip $@
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# further tools needing specialised compilation options
# -----------------------------------------------------------------------------
bin/heliphase: src/calc/heliphase.cpp src/core/heli.cpp src/core/magsys.cpp ext/tlibs2/libs/log.cpp
	$(CXX) $(STD) $(OPT) $(INCS) -DDEF_HELI_ORDER=4 -DNO_REDEFINITIONS -D__HACK_FULL_INST__ $(LIBDIRS) -o $@ $+ -lMinuit2 -llapacke -lgomp
	strip $@

bin/skx_gs: src/calc/skx_gs.cpp src/core/skx.cpp src/core/heli.cpp src/core/magsys.cpp ext/tlibs2/libs/log.cpp
	$(CXX) $(STD) $(OPT) $(INCS) -DDEF_SKX_ORDER=7 -DDEF_HELI_ORDER=7 -DNO_REDEFINITIONS -D__HACK_FULL_INST__ $(LIBDIRS) -o $@ $+ -lMinuit2 -llapacke -lgomp
	strip $@
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# general rules
# -----------------------------------------------------------------------------
%.o: %.cpp
	@echo "Compiling $< -> $@..."
	$(CXX) $(STD) $(OPT) $(DEFS) $(INCS) $(LIBDEFS) -c $< -o $@
# -----------------------------------------------------------------------------
