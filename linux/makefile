# Makefile for Linux mprime
#
# NOTE: To convert the included object files to ELF, your binutils 
# installation must be able to read the coff-i386 object format.  
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

FACTOROBJ = factor32.o
DUMMYOBJS = dummy0.o dummy4.o dummy8.o dummy12.o dummy16.o dummy20.o dummy24.o dummy28.o
LINUXOBJS = prime.o menu.o
EXE      = mprime
EXE2	 = sprime

# Insert proper dummy file to get _GWDATA data segment aligned
# on a 32 byte boundary.  Use proper dummyNN.o to do the alignment
# This is no longer necessary thanks to Alex Kruppa's gwnum.ld trick!

DUMMY1   =
DUMMY2   =

#########################################################################

$(EXE): security $(LINUXOBJS) $(FACTOROBJ) $(DUMMYOBJS)
	$(CC) $(LFLAGS) -o $(EXE) $(LINUXOBJS) $(DUMMY1) $(FACTOROBJ) $(LIBS)

$(EXE2): security $(LINUXOBJS) $(FACTOROBJ) $(DUMMYOBJS)
	$(CC) $(LFLAGS) -o $(EXE2) $(LINUXOBJS) $(DUMMY2) $(FACTOROBJ) $(LIBS2)

security:
	[ ! -e ../security.h ] && touch ../security.h || true
	[ ! -e ../security.c ] && touch ../security.c || true

clean:
	rm -f $(EXE) $(EXE2) $(LINUXOBJS) $(FACTOROBJ) $(DUMMYOBJS)

$(FACTOROBJ) $(DUMMYOBJS): %.o : %.obj
	objcopy -Oelf32-i386 -Icoff-i386 $< $@

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
