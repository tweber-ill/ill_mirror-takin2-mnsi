#
# MnSi dynamics module for Takin and helper tools
# @author Tobias Weber <tweber@ill.fr>
# @date 2018-2020
# @license GPLv2 (see 'LICENSE' file)
#

mingw_build = 0
debug_build = 0
strip_bins = 1


# -----------------------------------------------------------------------------
# setup
# -----------------------------------------------------------------------------
ifneq ($(mingw_build), 1)
	ifeq ("$(CXX)", "")
		CXX = g++
	endif

	SYSINCS = -I/usr/local/include \
		-I/usr/include/lapacke -I/usr/local/opt/lapack/include \
		-I/usr/include/qt5 -I/usr/include/x86_64-linux-gnu/qt5/ \
		-I/usr/local/include/Minuit2 \
		#-I/usr/local/Cellar/qt/5.15.0/include \
		#-I/home/tw/build/boost_1_73_0
	LIBDIRS = -L/usr/local/opt/lapack/lib -L/usr/local/lib

	LIBBOOSTSYS = -lboost_system
	LIBBOOSTFILESYS = -lboost_filesystem
	LIBBOOSTIO = -lboost_iostreams
	LIBBOOSOPTS = -lboost_program_options

	BIN_SUFFIX =
else
	CXX = x86_64-w64-mingw32-g++

	SYSINCS = -I/usr/x86_64-w64-mingw32/sys-root/mingw/include \
		-I/usr/x86_64-w64-mingw32/sys-root/mingw/include/qt5 \
		-I/usr/x86_64-w64-mingw32/sys-root/mingw/include/Minuit2
	LIBDIRS = -L/usr/x86_64-w64-mingw32/sys-root/mingw/bin/

	LIBBOOSTSYS = -lboost_system-x64
	LIBBOOSTFILESYS = -lboost_filesystem-x64
	LIBBOOSTIO = -lboost_iostreams-x64
	LIBBOOSOPTS = -lboost_program_options-x64

	BIN_SUFFIX = .exe
endif


ifneq ($(debug_build), 1)
	OPT = -O2 -w #-march=native

	ifeq ($(strip_bins), 1)
		STRIP = strip
	else
		STRIP = echo -e "Not stripping"
	endif
else
	OPT = -g -ggdb -Wall -Wextra
	STRIP = echo -e "Not stripping"
endif


#
# Setting DEF_SKX_ORDER to a higher value makes the backfolded spectrum for
# q perpendicular to B more accurate for high qs beyond the first magnetic
# Brillouin zone. The loss of precision at lower orders can be seen with the
# magnon energy levels not completely connecting at the zone boundaries for
# large qs and Es.
#
STD = -std=c++17
LIBDEFS = -fPIC
DEFS = -DDEF_SKX_ORDER=7 -DDEF_HELI_ORDER=7 \
	-DNO_MINIMISATION -DNO_REDEFINITIONS \
	-DSKX_USE_HOC=0 -DHELI_USE_HOC=1 \
	-D__HACK_FULL_INST__ #-DPLUGIN_APPLI
INCS = -Isrc -Iext -Iext/takin $(SYSINCS)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# meta rules
# -----------------------------------------------------------------------------
.PHONY: all clean

all: prepare \
	bin/genskx bin/genheli bin/merge bin/convert bin/dump \
	bin/drawskx bin/dyn bin/weight \
	bin/tof_img bin/tof_mask bin/tof_pol bin/tof_unite \
	bin/heliphase bin/heli_gs bin/skx_gs bin/weight_sum \
	lib/skxmod.so lib/skxmod_grid.so

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
lib/skxmod.so: src/takin/takin.o src/core/skx.o src/core/fp.o src/core/heli.o \
		src/core/longfluct.o src/core/magsys.o \
		ext/takin/tools/monteconvo/sqwbase.o \
		ext/tlibs2/libs/log.o ext/tlibs/log/log.o
	@echo "Linking Takin module $@..."
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -shared -o $@ $+ \
		-llapacke
	$(STRIP) $@

lib/skxmod_grid.so: src/takin/takin_grid.o \
		ext/takin/tools/monteconvo/sqwbase.o \
		ext/tlibs2/libs/log.o ext/tlibs/log/log.o
	@echo "Linking Takin grid module $@..."
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -shared -o $@ $+ \
		$(LIBBOOSTSYS) -lQt5Core
	$(STRIP) $@
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# tools
# -----------------------------------------------------------------------------
bin/genskx: src/takin/genskx.o src/core/skx.o \
		src/core/magsys.o ext/tlibs2/libs/log.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -o $@ $+ $(LIBBOOSTFILESYS) -llapacke -lpthread
	$(STRIP) $@$(BIN_SUFFIX)

bin/genheli: src/takin/genheli.o src/core/heli.o \
		src/core/magsys.o ext/tlibs2/libs/log.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -o $@ $+ $(LIBBOOSTFILESYS) -llapacke -lpthread
	$(STRIP) $@$(BIN_SUFFIX)

bin/merge: src/takin/merge.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+ $(LIBBOOSTSYS)
	$(STRIP) $@$(BIN_SUFFIX)

bin/convert: src/takin/convert.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -o $@ $+ $(LIBBOOSOPTS)
	$(STRIP) $@$(BIN_SUFFIX)

bin/dump: src/takin/dump.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+
	$(STRIP) $@$(BIN_SUFFIX)

bin/drawskx: src/calc/drawskx.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+
	$(STRIP) $@$(BIN_SUFFIX)

bin/tof_pol: src/calc/tof_pol.o ext/tlibs2/libs/log.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+ \
		$(LIBBOOSTSYS) $(LIBBOOSTFILESYS) $(LIBBOOSTIO) -lpng -lMinuit2 -lMinuit2Math
	$(STRIP) $@$(BIN_SUFFIX)

bin/tof_unite: src/calc/tof_unite.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+ \
		$(LIBBOOSTSYS) $(LIBBOOSTFILESYS) $(LIBBOOSTIO)
	$(STRIP) $@$(BIN_SUFFIX)

bin/tof_mask: src/calc/tof_mask.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+ \
		$(LIBBOOSTSYS) $(LIBBOOSTFILESYS) $(LIBBOOSTIO) -lpng
	$(STRIP) $@$(BIN_SUFFIX)

bin/tof_img: src/calc/tof_img.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) -o $@ $+ \
		$(LIBBOOSTSYS) $(LIBBOOSTFILESYS) $(LIBBOOSTIO) -lpng
	$(STRIP) $@$(BIN_SUFFIX)

bin/dyn: src/calc/dyn.o src/core/skx.o src/core/fp.o src/core/heli.o \
		src/core/magsys.o ext/tlibs2/libs/log.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -o $@ $+ $(LIBBOOSOPTS) -llapacke -lpthread
	$(STRIP) $@$(BIN_SUFFIX)

bin/weight: src/calc/weight.o src/core/skx.o src/core/fp.o src/core/heli.o \
		src/core/magsys.o ext/tlibs2/libs/log.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -o $@ $+ -llapacke -lpthread
	$(STRIP) $@$(BIN_SUFFIX)

bin/weight_sum: src/calc/weight_sum.o src/core/skx.o src/core/fp.o src/core/heli.o \
		src/core/magsys.o ext/tlibs2/libs/log.o
	$(CXX) $(STD) $(OPT) $(DEFS) $(LIBDIRS) $(LIBDEFS) -o $@ $+ -llapacke -lpthread
	$(STRIP) $@$(BIN_SUFFIX)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# further tools needing specialised compilation options
# -----------------------------------------------------------------------------
bin/heliphase: src/calc/heliphase.cpp src/core/heli.cpp \
		src/core/magsys.cpp ext/tlibs2/libs/log.cpp
	$(CXX) $(STD) $(OPT) $(INCS) -DDEF_HELI_ORDER=4 \
		-DNO_REDEFINITIONS -D__HACK_FULL_INST__ \
		$(LIBDIRS) -o $@ $+ -lMinuit2 -lMinuit2Math -llapacke
	$(STRIP) $@$(BIN_SUFFIX)

bin/heli_gs: src/calc/heli_gs.cpp src/core/heli.cpp \
		src/core/magsys.cpp ext/tlibs2/libs/log.cpp
	$(CXX) $(STD) $(OPT) $(INCS) -DDEF_HELI_ORDER=9 \
		-DNO_REDEFINITIONS -D__HACK_FULL_INST__ \
		$(LIBDIRS) -o $@ $+ -lMinuit2 -lMinuit2Math -llapacke
	$(STRIP) $@$(BIN_SUFFIX)

bin/skx_gs: src/calc/skx_gs.cpp src/core/skx.cpp src/core/heli.cpp \
		src/core/magsys.cpp ext/tlibs2/libs/log.cpp
	$(CXX) $(STD) $(OPT) $(INCS) -DDEF_SKX_ORDER=9 -DDEF_HELI_ORDER=9 \
		-DNO_REDEFINITIONS -D__HACK_FULL_INST__ \
		$(LIBDIRS) -o $@ $+ -lMinuit2 -lMinuit2Math -llapacke
	$(STRIP) $@$(BIN_SUFFIX)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# general rules
# -----------------------------------------------------------------------------
%.o: %.cpp
	@echo "Compiling $< -> $@..."
	$(CXX) $(STD) $(OPT) $(DEFS) $(INCS) $(LIBDEFS) -c $< -o $@
# -----------------------------------------------------------------------------
