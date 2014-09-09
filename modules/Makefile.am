## Process this file with automake to produce Makefile.in

SUBDIRS = RStarTree kdt fes2004 wavewatch culvert

AM_CPPFLAGS = -DGFS_MODULES_DIR=\"$(libdir)/gerris\"

INCLUDES = -I$(top_srcdir)/src -I$(includedir) \
           -DG_LOG_DOMAIN=\"Gfs-modules\" $(GTS_CFLAGS)

if HAS_LIBPROJ
MAP = libmap2D.la libmap3D.la
endif
if BUILD_TIDE
TIDE = libtide2D.la libtide3D.la
endif
if BUILD_STOKES
STOKES = libstokes2D.la libstokes3D.la
endif
if BUILD_WAVEWATCH
WAVEWATCH = libwavewatch2D.la
endif
if BUILD_HYPRE
HYPRE = libhypre2D.la libhypre3D.la
endif
if BUILD_LIS
LIS = liblis2D.la liblis3D.la
endif
if BUILD_AGMG
AGMG = libagmg2D.la libagmg3D.la
endif
if BUILD_ODE
ODE = libode2D.la libode3D.la
endif
if HAS_GSL
BUBBLES = libbubbles2D.la libbubbles3D.la
endif
if HAS_FFTW3
FFT = libfft2D.la libfft3D.la
endif

pkglib_LTLIBRARIES = \
	$(MAP) \
	libterrain2D.la \
	libterrain3D.la \
	libparticulates2D.la \
	libparticulates3D.la \
	libokada2D.la \
	libokada3D.la \
	libdf33D.la \
	libelectrohydro2D.la \
	libelectrohydro3D.la \
	libskewsymmetric2D.la \
	libskewsymmetric3D.la \
	libtopics2D.la \
	libtopics3D.la \
	libculvert2D.la \
	liblayered2D.la \
	$(TIDE) \
	$(STOKES) \
	$(WAVEWATCH) \
	$(HYPRE) \
	$(LIS) \
	$(AGMG) \
	$(ODE) \
        $(BUBBLES) \
	$(FFT)

bin_PROGRAMS = \
	rsurface2kdt

AM_LDFLAGS = $(NO_UNDEFINED)\
        -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic -module

libmap3D_la_SOURCES = map.c
libmap3D_la_LIBADD = $(GFS3D_LIBS) -lproj
libmap2D_la_SOURCES = map.c
libmap2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libmap2D_la_LIBADD = $(GFS2D_LIBS) -lproj

KDTLIBS = -Lkdt -lkdt -lm
KDTDEPS = kdt/libkdt.la

libterrain3D_la_SOURCES = terrain.c
libterrain3D_la_LIBADD = $(GFS3D_LIBS) $(KDTLIBS)
libterrain3D_la_DEPENDENCIES = $(KDTDEPS)

libterrain2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libterrain2D_la_SOURCES = terrain.c
libterrain2D_la_LIBADD = $(GFS2D_LIBS) $(KDTLIBS)
libterrain2D_la_DEPENDENCIES = $(KDTDEPS)

libtide3D_la_SOURCES = tide.c
libtide3D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS)
libtide3D_la_LIBADD = $(GFS3D_LIBS) -Lfes2004 -lfes -lnetcdf $(GSL_LIBS)
libtide2D_la_SOURCES = tide.c
libtide2D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) -DFTT_2D=1
libtide2D_la_LIBADD = $(GFS2D_LIBS) -Lfes2004 -lfes -lnetcdf $(GSL_LIBS)

libstokes3D_la_SOURCES = stokes.c CW263.f
libstokes3D_la_CFLAGS = $(AM_CFLAGS)
libstokes3D_la_LIBADD = $(GFS3D_LIBS)
libstokes2D_la_SOURCES = stokes.c CW263.f
libstokes2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libstokes2D_la_LIBADD = $(GFS2D_LIBS)

libwavewatch2D_la_SOURCES = wavewatch.c
libwavewatch2D_la_LIBADD = $(GFS2D_LIBS) -Lwavewatch -lwavewatch -L/usr/lib/gcc/i486-linux-gnu/4.2 -lgfortran
libwavewatch2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1 \
  -DW3SRCE=`nm wavewatch/libwavewatch.a | grep w3srce | awk '{ if ($$2 == "T") print $$3; }'` \
  -DGFSW3INIT=`nm wavewatch/libwavewatch.a | grep gfsw3_init | awk '{ if ($$2 == "T") print $$3; }'`
libwavewatch2D_la_DEPENDENCIES = wavewatch/libwavewatch.a

libode3D_la_SOURCES = ode.c
libode3D_la_CFLAGS = $(AM_CFLAGS) $(ODE_CFLAGS)
libode3D_la_LIBADD = $(GFS3D_LIBS) $(ODE_LIBS)
libode2D_la_SOURCES = ode.c
libode2D_la_CFLAGS = $(AM_CFLAGS) $(ODE_CFLAGS) -DFTT_2D=1
libode2D_la_LIBADD = $(GFS2D_LIBS) $(ODE_LIBS)

rsurface2kdt_SOURCES = rsurface2kdt.c rsurface.c rsurface.h
rsurface2kdt_LDADD = -LRStarTree -lcSmRST $(KDTLIBS)
rsurface2kdt_CFLAGS = $(AM_CFLAGS)
rsurface2kdt_DEPENDENCIES = $(KDTDEPS)

noinst_LTLIBRARIES = libparticulatecommon2D.la libparticulatecommon3D.la

libparticulatecommon3D_la_SOURCES = particulatecommon.c particulatecommon.h
libparticulatecommon3D_la_CFLAGS = $(AM_CFLAGS)
libparticulatecommon2D_la_SOURCES = particulatecommon.c particulatecommon.h
libparticulatecommon2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1

libparticulates3D_la_SOURCES = particulates.c
libparticulates3D_la_CFLAGS = $(AM_CFLAGS)
libparticulates3D_la_LIBADD = libparticulatecommon3D.la $(GFS3D_LIBS)
libparticulates3D_la_DEPENDENCIES = libparticulatecommon3D.la
libparticulates2D_la_SOURCES = particulates.c
libparticulates2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libparticulates2D_la_LIBADD = libparticulatecommon2D.la $(GFS2D_LIBS)
libparticulates2D_la_DEPENDENCIES = libparticulatecommon2D.la

libbubbles3D_la_SOURCES = bubbles.c
libbubbles3D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS)
libbubbles3D_la_LIBADD = libparticulatecommon3D.la $(GFS3D_LIBS) $(GSL_LIBS)
libbubbles3D_la_DEPENDENCIES = libparticulatecommon3D.la
libbubbles2D_la_SOURCES = bubbles.c
libbubbles2D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) -DFTT_2D=1
libbubbles2D_la_LIBADD = libparticulatecommon2D.la $(GFS2D_LIBS) $(GSL_LIBS)
libbubbles2D_la_DEPENDENCIES = libparticulatecommon2D.la

libfft3D_la_SOURCES = fft.c
libfft3D_la_CFLAGS = $(AM_CFLAGS) $(FFTW3_CFLAGS)
libfft3D_la_LIBADD = $(GFS3D_LIBS) $(FFTW3_LIBS)
libfft2D_la_SOURCES = fft.c
libfft2D_la_CFLAGS = $(AM_CFLAGS) $(FFTW3_CFLAGS) -DFTT_2D=1
libfft2D_la_LIBADD = $(GFS2D_LIBS) $(FFTW3_LIBS)

libhypre3D_la_SOURCES = hypre.c
libhypre3D_la_CFLAGS = $(AM_CFLAGS)
libhypre3D_la_LIBADD = $(GFS3D_LIBS) -lHYPRE -lm
libhypre2D_la_SOURCES = hypre.c
libhypre2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libhypre2D_la_LIBADD = $(GFS2D_LIBS) -lHYPRE -lm

liblis3D_la_SOURCES = lis.c
liblis3D_la_CFLAGS = $(AM_CFLAGS)
liblis3D_la_LIBADD = $(GFS3D_LIBS) -llis -lgfortran -lm
liblis2D_la_SOURCES = lis.c
liblis2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
liblis2D_la_LIBADD = $(GFS2D_LIBS) -llis -lgfortran -lm

libagmg3D_la_SOURCES = agmg.c dagmg.f90
libagmg3D_la_CFLAGS = $(AM_CFLAGS)
libagmg3D_la_LIBADD = $(GFS3D_LIBS) -llapack -lblas -lgfortran -lm
libagmg2D_la_SOURCES = agmg.c dagmg.f90
libagmg2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libagmg2D_la_LIBADD = $(GFS2D_LIBS) -llapack -lblas -lgfortran -lm

libokada3D_la_SOURCES = okada.c
libokada3D_la_CFLAGS = $(AM_CFLAGS)
libokada3D_la_LIBADD = $(GFS3D_LIBS)
libokada2D_la_SOURCES = okada.c
libokada2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libokada2D_la_LIBADD = $(GFS2D_LIBS)

libdf33D_la_SOURCES = df3.c
libdf33D_la_CFLAGS = $(AM_CFLAGS)
libdf33D_la_LIBADD = $(GFS3D_LIBS)

libelectrohydro3D_la_SOURCES = electrohydro.c
libelectrohydro3D_la_CFLAGS = $(AM_CFLAGS)
libelectrohydro3D_la_LIBADD = $(GFS3D_LIBS)
libelectrohydro2D_la_SOURCES = electrohydro.c
libelectrohydro2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libelectrohydro2D_la_LIBADD = $(GFS2D_LIBS)

libskewsymmetric3D_la_SOURCES = skewsymmetric.c
libskewsymmetric3D_la_CFLAGS = $(AM_CFLAGS)
libskewsymmetric3D_la_LIBADD = $(GFS3D_LIBS)
libskewsymmetric2D_la_SOURCES = skewsymmetric.c
libskewsymmetric2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libskewsymmetric2D_la_LIBADD = $(GFS2D_LIBS)

libtopics3D_la_SOURCES = topics.c
libtopics3D_la_CFLAGS = $(AM_CFLAGS)
libtopics3D_la_LIBADD = $(GFS3D_LIBS)
libtopics2D_la_SOURCES = topics.c
libtopics2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libtopics2D_la_LIBADD = $(GFS2D_LIBS)

libculvert2D_la_SOURCES = culvert.c
libculvert2D_la_CFLAGS = $(AM_CFLAGS)
libculvert2D_la_LIBADD = $(GFS2D_LIBS) -Lculvert -lboyd87 -lm
libculvert2D_la_DEPENDENCIES = culvert/libboyd87.la

liblayered2D_la_SOURCES = layered.c
liblayered2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
liblayered2D_la_LIBADD = $(GFS2D_LIBS)
