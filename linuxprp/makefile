# Makefile for Linux PRP
#
# NOTE: To convert the included object files to ELF, your binutils 
# installation must be able to read the coff-i386 object format.  
#
# Assuming your binutils are dynamically linked, you can add COFF support
# by configuring binutils with the switch '--enable-targets=i386-go32'
# and rebuilding libbfd.so .

CC = gcc
CFLAGS = -I.. -O2 -m486 -malign-double

#CFLAGS for Pentium (I don't know who recommended these settings)
#CFLAGS = -I.. -O2 -malign-double -malign-loops=2 -malign-jumps=2 -malign-functions=2

# Setting suggested by FenwayFrank - I don't know what it does
VPATH=../prime95

LFLAGS = -Wl,-M
LIBS   = -lm
LIBS2  = -lm -static

COFFOBJS = ../linux/mult.o ../linux/mult1.o ../linux/mult2.o ../linux/mult2a.o ../linux/mult3.o ../linux/mult3a.o ../linux/mult4.o ../linux/mult4a.o ../linux/mult4b.o ../linux/mult1p.o ../linux/mult2p.o ../linux/mult2ap.o ../linux/mult3p.o ../linux/mult3ap.o ../linux/mult3q.o ../linux/mult3aq.o ../linux/mult4p.o ../linux/mult4ap.o ../linux/mult4bp.o ../linux/mult4q.o ../linux/mult4aq.o ../linux/mult4bq.o ../linux/mult1aux.o ../linux/mult2aux.o ../linux/mult3aux.o ../linux/mult3auq.o ../linux/mult4aux.o ../linux/mult4auq.o ../linux/xmult1.o ../linux/xmult1ax.o ../linux/xmult2.o ../linux/xmult2a.o ../linux/xmult2ax.o ../linux/xmult3.o ../linux/xmult3a.o ../linux/xmult3ax.o ../linux/ecmhelp.o ../linux/cpuid.o
DUMMYOBJS = ../linux/dummy4.o ../linux/dummy8.o ../linux/dummy12.o ../linux/dummy16.o ../linux/dummy20.o ../linux/dummy24.o ../linux/dummy28.o
LINUXOBJS = prime.o menu.o
EXE      = prp
EXE2	 = sprp
DUMMY1   = ../linux/dummy24.o
DUMMY2   = ../linux/dummy24.o

#########################################################################

$(EXE): $(LINUXOBJS) $(COFFOBJS) $(DUMMYOBJS)
	$(CC) $(LFLAGS) -o $(EXE) $(LINUXOBJS) $(DUMMY1) $(COFFOBJS) $(LIBS)

$(EXE2): $(LINUXOBJS) $(COFFOBJS) $(DUMMYOBJS)
	$(CC) $(LFLAGS) -o $(EXE2) $(LINUXOBJS) $(DUMMY2) $(COFFOBJS) $(LIBS2)

clean:
	rm -f $(EXE) $(EXE2) $(LINUXOBJS) $(COFFOBJS) $(DUMMYOBJS)

.c.o:
	$(CC) $(CFLAGS) -c $<
