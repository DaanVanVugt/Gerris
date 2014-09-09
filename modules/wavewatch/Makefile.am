## Process this file with automake to produce Makefile.in

if BUILD_WAVEWATCH
WAVEWATCH = libwavewatch.a
CHECKTYPES = checktypes
endif

noinst_LIBRARIES = $(WAVEWATCH)
noinst_PROGRAMS = $(CHECKTYPES)

EXTRA_DIST = gfsad3 \
	w3init222.ftn wavewatch_222.h \
	wavewatch_312.h \
	w3init314.ftn wavewatch_314.h \
	fchecktypes.ftn

checktypes_SOURCES = cchecktypes.c
checktypes_LDADD = fchecktypes.o

libwavewatch_a_SOURCES =

CLEANFILES = cfortrantypes.h

libwavewatch.a: w3initmd.o
	rm -f libwavewatch.a
	$(AR) cru libwavewatch.a `grep WWATCH3_DIR $(HOME)/.wwatch3.env | awk '{print $$2 "/obj/w3*.o"}'` w3initmd.o
	$(RANLIB) libwavewatch.a

w3initmd.o: $(W3INIT).ftn $(HOME)/.wwatch3.env gfsad3 cfortrantypes.h
	sh gfsad3 $(W3INIT).ftn
	mv -f $(W3INIT).o w3initmd.o

fchecktypes.o: fchecktypes.ftn $(HOME)/.wwatch3.env gfsad3
	sh gfsad3 fchecktypes.ftn

cfortrantypes.h: checktypes
	./checktypes > cfortrantypes.h

w3init312.ftn: w3init314.ftn
	cp -f w3init314.ftn w3init312.ftn
