# Makefile for Linux gwnum library
#
# NOTE: To convert the included object files to ELF, your binutils 
# installation must be able to read the coff-i386 object format.  
#
# Assuming your binutils are dynamically linked, you can add COFF support
# by configuring binutils with the switch '--enable-targets=i386-go32'
# and rebuilding libbfd.so .

CC = gcc
CFLAGS = -I.. -O2 -march=i486 -malign-double -DADD_UNDERSCORES

CPP = g++
CPPFLAGS = -I.. -I../qd -O2 -march=i486 -malign-double -DADD_UNDERSCORES

COFFOBJS = mult.o mult1.o mult2.o mult2a.o mult3.o mult3a.o mult4.o mult2p.o mult2ap.o mult3p.o mult3ap.o mult4p.o mult1aux.o mult2aux.o xmult1.o xmult1ax.o xmult2.o xmult2a.o xmult2ax.o xmult3.o xmult3a.o xmult4.o xmult5.o xmult6.o xmult2p.o xmult2ap.o xmult3p.o xmult3ap.o xmult4p.o xmult5p.o xmult6p.o gianthlp.o cpuidhlp.o timeit.o timeitp.o
LINUXOBJS = cpuid.o gwnum.o gwutil.o gwdbldbl.o giants.o ecmstag1.o

LIB = gwnum.lib

#########################################################################

$(LIB): $(LINUXOBJS) $(COFFOBJS)
	rm -f gwnum.lib
	ar -qs gwnum.lib $(LINUXOBJS) $(COFFOBJS)

clean:
	rm -f $(LIB) $(LINUXOBJS) $(COFFOBJS)

$(COFFOBJS): %.o : %.obj
	objcopy -Oelf32-i386 -Icoff-i386 $< $@

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
