# Makefile for Mac OS X 64-bit mprime
#
# Hwloc can be downloaded from https://www.open-mpi.org/projects/hwloc/
# Configure hwloc with --enable-static --disable-shared --disable-io
#
# GMP can be configured with --prefix=/usr/local --enable-fat --enable-static --disable-shared
#
# If troubles occur loading <atomic>
# Try installing gcc-12 using homebrew, then setting
# CC = /usr/local/bin/gcc-12
# CPP = /usr/local/bin/g++-12
# add -I/usr/local/include/ to CFLAGS and CPPFLAGS 

ENVP= MACOSX_DEPLOYMENT_TARGET=10.9

CC = gcc
CFLAGS = -I.. -I../gwnum -I../linux -std=gnu99 -O2 -DX86_64 -DCOMMAND_LINE_MPRIME -m64 -DPORT=10

CPP = g++
CPPFLAGS = -I.. -I../gwnum -I../linux -I../../boost_1_76_0 -std=c++11 -O2 -DX86_64 -DCOMMAND_LINE_MPRIME -m64 -DPORT=10

LFLAGS = -m64 -Wl,-no_pie -L/usr/local/lib
LIBS   = ../gwnum/gwnum.a ../gwnum/polymult.a -lm -lpthread /usr/local/lib/libhwloc.a /usr/local/lib/libgmp.a -lcurl -framework IOKit -framework CoreFoundation -lc++

FACTOROBJ = ../prime95/macosx64/factor64.o
OBJS = prime.o menu.o cJSON.o ecm.o exponentiate.o pair.o pm1prob.o
EXE = mprime

#########################################################################

$(EXE): security $(OBJS) $(FACTOROBJ)
	$(ENVP) $(CC) $(LFLAGS) -o $(EXE) $(OBJS) $(FACTOROBJ) $(LIBS)

security:
	[ ! -e ../security.h ] && touch ../security.h || true
	[ ! -e ../security.c ] && touch ../security.c || true
	[ ! -e ../secure5.c ] && touch ../secure5.c || true

clean:
	rm -f $(EXE) $(EXE2) $(OBJS)

cJSON.o:
	$(ENVP) $(CC) $(CFLAGS) -c ../cJSON.c

ecm.o:
	$(ENVP) $(CPP) $(CPPFLAGS) -c ../ecm.cpp

exponentiate.o:
	$(ENVP) $(CC) $(CFLAGS) -c ../exponentiate.c

pair.o:
	$(ENVP) $(CPP) $(CPPFLAGS) -c ../pair.cpp

pm1prob.o:
	$(ENVP) $(CC) $(CFLAGS) -c ../pm1prob.c

.c.o:
	$(ENVP) $(CC) $(CFLAGS) -c $<

.cpp.o:
	$(ENVP) $(CPP) $(CPPFLAGS) -c $<
