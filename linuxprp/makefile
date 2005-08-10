# Makefile for Linux PRP
#
# Assuming your binutils are dynamically linked, you can add COFF support
# by configuring binutils with the switch '--enable-targets=i386-go32'
# and rebuilding libbfd.so .

CC = gcc
CFLAGS = -I.. -I../gwnum -O2 -march=i486 -malign-double

CPP = g++
CPPFLAGS = -I.. -I../gwnum -O2 -march=i486 -malign-double

LFLAGS = -Wl,-M
LIBS   = ../gwnum/gwnum.lib ../gwnum/gwnum.ld -lm -lstdc++
LIBS2  = ../gwnum/gwnum.lib ../gwnum/gwnum.ld -lm -lstdc++ -static

LINUXOBJS = prime.o menu.o
EXE      = prp
EXE2	 = sprp

# Insert proper dummy file to get _GWDATA data segment aligned
# on a 32 byte boundary.  Use proper dummyNN.o to do the alignment
# This is no longer necessary thanks to Alex Kruppa's gwnum.ld trick!

DUMMY1   =
DUMMY2   =

#########################################################################

$(EXE): $(LINUXOBJS)
	$(CC) $(LFLAGS) -o $(EXE) $(LINUXOBJS) $(DUMMY1) $(LIBS)

$(EXE2): $(LINUXOBJS)
	$(CC) $(LFLAGS) -o $(EXE2) $(LINUXOBJS) $(DUMMY2) $(LIBS2)

clean:
	rm -f $(EXE) $(EXE2) $(LINUXOBJS)

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
