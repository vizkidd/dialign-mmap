###############################
#                             #
#   Makefile for dialign-mmap #
#                             #
###############################
#
#  make            build ./dialign2-2  (OpenMP + mmap store)
#  make NOOMP=1    build without OpenMP
#  make DEBUG=1    build with -O0 -g and AddressSanitizer/UBSan
#
CC      ?= gcc
# -ffp-contract=off : no fused multiply-add, so that float results are
#                     bit-identical on every platform / compiler setting
OPT     ?= -O2
CFLAGS  = -c $(OPT) -g -std=gnu99 -ffp-contract=off -DCONS -I.
LIBS    = -lm -lpthread

ifndef NOOMP
CFLAGS += -fopenmp
LIBS   += -lgomp
endif

ifdef DEBUG
CFLAGS  = -c -O0 -g -std=gnu99 -ffp-contract=off -DCONS -I. -fsanitize=address,undefined -fno-omit-frame-pointer
LIBS    = -lm -lpthread -fsanitize=address,undefined
endif

OBJS = dialign.o functions.o para.o \
       output.o wgt.o anchor.o regex.o \
       input.o frag_chain.o \
       alig_graph_closure.o pratique.o mmstore.o frags.o

dialign2-2: $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o dialign2-2

%.o: %.c dialign.h define.h mmstore.h pratique.h alig_graph_closure.h
	$(CC) $(CFLAGS) $<

clean:
	rm -f *.o dialign2-2
