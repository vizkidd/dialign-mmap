########################################################################
#  Top-level Makefile for dialign-mmap
#
#  make            build ./src/dialign2-2 with auto-detected build flags
#  make install    install the binary, a DIALIGN2_DIR-setting wrapper,
#                  and the reference data files (PREFIX=/usr/local default)
#  make test       run the self-contained test suite (no reference
#                  original needed)
#  make test ORIG=/path/to/original/dialign2-2 \
#            ORIG_SRC=/path/to/original/src \
#            run the full differential suite as well (matrix,
#            robustness, determinism, tree unit test, fuzz)
#  make clean / distclean / uninstall
#
#  Override variables on the command line, e.g.:
#    make PREFIX=$HOME/.local install
#    make CC=clang
#    make NOOMP=1              force-disable OpenMP
#    make OPT=-O2              force a specific optimization level
########################################################################

SHELL   := /bin/bash
ROOT    := $(CURDIR)
SRC_DIR := $(ROOT)/src
TESTS   := $(ROOT)/tests
DATA_DIR:= $(ROOT)/dialign2_dir
BIN     := dialign2-2

ifeq ($(origin CC),default)
  CC := gcc
endif

PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin
DATADIR ?= $(PREFIX)/share/dialign-mmap
DESTDIR ?=

# ----------------------------------------------------------------------
# Auto-detect OpenMP support and whether -O3 is safe on this compiler /
# codebase (verified byte-identical to -O2, not just "the flag is
# accepted" -- see tools/detect.sh). Both NOOMP=1 and OPT=... on the
# command line override the corresponding auto-detected value; the
# detection only runs (and only costs a couple of test-alignment runs)
# when the corresponding variable wasn't already set.
# ----------------------------------------------------------------------
ifeq ($(origin NOOMP),undefined)
  ifeq ($(origin OPT),undefined)
    DETECTED       := $(shell $(ROOT)/tools/detect.sh $(CC) $(ROOT) 2>/dev/null)
    DETECTED_OPENMP:= $(filter OPENMP=1,$(DETECTED))
    OPT            := $(patsubst OPT=%,%,$(filter OPT=%,$(DETECTED)))
  else
    DETECTED_OPENMP:= $(shell printf '\#include <omp.h>\nint main(){return omp_get_num_threads()>=0?0:1;}\n' > /tmp/.dmm_omp_$$$$.c; $(CC) -fopenmp /tmp/.dmm_omp_$$$$.c -o /tmp/.dmm_omp_$$$$ -lgomp >/dev/null 2>&1 && echo OPENMP=1; rm -f /tmp/.dmm_omp_$$$$.c /tmp/.dmm_omp_$$$$)
  endif
  ifeq ($(DETECTED_OPENMP),)
    NOOMP := 1
  endif
endif
OPT ?= -O2

.PHONY: all build install uninstall test test-full clean distclean print-config

all: build

build:
	@echo "== dialign-mmap: CC=$(CC)  OPT=$(OPT)  OpenMP=$(if $(NOOMP),no,yes) =="
	$(MAKE) -C $(SRC_DIR) CC="$(CC)" OPT="$(OPT)" $(if $(NOOMP),NOOMP=1)
	@echo "== built $(SRC_DIR)/$(BIN) =="

print-config:
	@echo "CC      = $(CC)"
	@echo "OPT     = $(OPT)"
	@echo "NOOMP   = $(NOOMP)"
	@echo "PREFIX  = $(PREFIX)"
	@echo "BINDIR  = $(DESTDIR)$(BINDIR)"
	@echo "DATADIR = $(DESTDIR)$(DATADIR)"

# ----------------------------------------------------------------------
# install: the real binary goes in as dialign2-2.bin; a small wrapper
# script named dialign2-2 is generated to export DIALIGN2_DIR (pointing
# at the installed reference-data copy) and exec the real binary, so
# `dialign2-2 seqs.fa` works right after install with no manual
# environment-variable setup (the single most common first-run error
# with the plain original tool).
# ----------------------------------------------------------------------
install: build
	install -d "$(DESTDIR)$(BINDIR)"
	install -d "$(DESTDIR)$(DATADIR)"
	install -m 755 "$(SRC_DIR)/$(BIN)" "$(DESTDIR)$(BINDIR)/$(BIN).bin"
	cp -R "$(DATA_DIR)/." "$(DESTDIR)$(DATADIR)/"
	printf '#!/bin/sh\nexec env DIALIGN2_DIR="%s" "%s/%s.bin" "$$@"\n' \
	  "$(DATADIR)" "$(BINDIR)" "$(BIN)" > "$(DESTDIR)$(BINDIR)/$(BIN)"
	chmod 755 "$(DESTDIR)$(BINDIR)/$(BIN)"
	@echo "== installed: $(DESTDIR)$(BINDIR)/$(BIN) (wrapper) =="
	@echo "==            $(DESTDIR)$(BINDIR)/$(BIN).bin (real binary) =="
	@echo "==            $(DESTDIR)$(DATADIR)/ (reference data, DIALIGN2_DIR) =="
	@echo "== run: $(BIN) <seqs.fa> =="

uninstall:
	rm -f  "$(DESTDIR)$(BINDIR)/$(BIN)" "$(DESTDIR)$(BINDIR)/$(BIN).bin"
	rm -rf "$(DESTDIR)$(DATADIR)"

# ----------------------------------------------------------------------
# test: always runs the self-contained checks, which need only the
# binary just built (new_options.sh: the port's added CLI options give
# identical output regardless of -nommap/-threads/-tmpdir; tree_unit.c:
# the rewritten tree-printing code agrees with a reference copy of the
# ORIGINAL av_tree_print(), so it needs ORIG_SRC but not a full original
# build).
#
# test-full (or `make test ORIG=...`) additionally runs everything that
# needs a full build of the plain original DIALIGN 2.2.1 as a reference:
# the option matrix, the 7-variant robustness sweep, determinism.sh, and
# a differential fuzz pass. Not part of this repository (the whole point
# of the port is to not require carrying a copy of the original around),
# so it only runs when you point ORIG (and, for the tree unit test,
# ORIG_SRC) at one.
# ----------------------------------------------------------------------
test: build
	@echo "== new_options.sh (self-contained) =="
	$(TESTS)/new_options.sh $(SRC_DIR)/$(BIN)
ifneq ($(ORIG_SRC),)
	@echo "== tree_unit.sh (needs ORIG_SRC's functions.c) =="
	$(TESTS)/tree_unit.sh $(ORIG_SRC)/functions.c
endif
ifneq ($(ORIG),)
	@$(MAKE) --no-print-directory test-full ORIG="$(ORIG)" ORIG_SRC="$(ORIG_SRC)"
else
	@echo
	@echo "== skipped matrix/robustness/determinism/fuzz: no ORIG= given =="
	@echo "== (make test ORIG=/path/to/original/dialign2-2 ORIG_SRC=/path/to/original/src for the full differential suite) =="
endif

test-full: build
	@test -n "$(ORIG)" || { echo "test-full needs ORIG=/path/to/original/dialign2-2"; exit 1; }
	@echo "== run_matrix.sh: $(words $(shell grep -vc '^\#' $(TESTS)/cases.txt)) cases vs. original =="
	$(TESTS)/run_matrix.sh "$(ORIG)" "$(SRC_DIR)/$(BIN)" "$(TESTS)/cases.txt"
	@echo "== determinism.sh =="
	$(TESTS)/determinism.sh "$(ORIG)" "$(SRC_DIR)/$(BIN)"
	@command -v python3 >/dev/null && { \
	  echo "== fuzz.py (100 runs, seed 1) =="; \
	  python3 $(TESTS)/fuzz.py "$(ORIG)" "$(SRC_DIR)/$(BIN)" 100 1; \
	} || echo "== skipped fuzz.py: python3 not found =="

clean:
	$(MAKE) -C $(SRC_DIR) clean

distclean: clean
	rm -f $(TESTS)/data/*.ali $(TESTS)/data/*.cap $(TESTS)/data/*.fop 2>/dev/null || true
