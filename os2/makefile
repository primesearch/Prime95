# Makefile for OS/2 mprime
# to compile with Watcom 11c / OpenWatcom
# use with GNUMAKE
# (c) Max Alekseyev

CC = wcc386
CFLAGS = /ox /d0 /I..;..\linux /5s /D__IBMC__ /DTCPV40HDRS
LINK = wlink
LFLAGS = LIBPath ..\prime95 op m import DosQuerySysState DOSCALLS.368

# debug version
# CFLAGS = /ox /d3 /I.. /5s /D__IBMC__ /D_DEBUG
# LFLAGS = LIBPath ..\prime95 op m debug all import DosQuerySysState DOSCALLS.368

COFFOBJS = cpuid.obj,ecmhelp.obj,factor64.obj,mult.obj,mult1.obj,mult1aux.obj,mult1p.obj,mult2.obj,mult2a.obj,mult2ap.obj,mult2aux.obj,mult2p.obj,mult3.obj,mult3a.obj,mult3ap.obj,mult3aq.obj,mult3auq.obj,mult3aux.obj,mult3p.obj,mult3q.obj,mult4.obj,mult4a.obj,mult4ap.obj,mult4aq.obj,mult4auq.obj,mult4aux.obj,mult4b.obj,mult4bp.obj,mult4bq.obj,mult4p.obj,mult4q.obj,xmult1.obj,xmult1ax.obj,xmult2.obj,xmult2a.obj,xmult2ax.obj,xmult3.obj,xmult3a.obj,xmult3ax.obj
DUMMYOBJS = dummy4.obj,dummy8.obj,dummy12.obj,dummy16.obj,dummy20.obj,dummy24.obj,dummy28.obj
PRIMEOBJS = prime.obj,menu.obj
# PRIMESRC = prime.c,menu.c
# LIBS = tcpip32.lib
LIBS = so32dll.lib,tcp32dll.lib
EXE = primeos2.exe

$(EXE): prime.obj menu.obj		# $(PRIMEOBJS),$(COFFOBJS),$(DUMMYOBJS)
	$(LINK) $(LFLAGS) file prime.obj,menu.obj name $(EXE) libfile $(COFFOBJS),dummy24.obj library $(LIBS) 

prime.obj: ..\linux\prime.c
	$(CC) $(CFLAGS) $<

menu.obj: ..\linux\menu.c
	$(CC) $(CFLAGS) $<

clean:
	del $(EXE) prime.obj menu.obj
