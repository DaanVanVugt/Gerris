## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = 1.6

AM_CPPFLAGS = -DGFS_MODULES_DIR=\"$(libdir)/gerris\" -DGFS_DATA_DIR=\"$(pkgdatadir)\"

INCLUDES = -I$(top_srcdir) -I$(includedir) \
           -DG_LOG_DOMAIN=\"Gfs\" $(GTS_CFLAGS)

gerris2D.pc: gerris2D.pc.in
gerris3D.pc: gerris3D.pc.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = gerris2D.pc gerris3D.pc
pkgdata_DATA = m4.awk

lib_LTLIBRARIES = libgfs2D.la libgfs3D.la

BUILT_SOURCES= \
	gerris2D.pc \
	gerris3D.pc

GFS_HDS = \
	ftt.h \
	fluid.h \
	variable.h \
	output.h \
	solid.h \
	poisson.h \
	advection.h \
	boundary.h \
	mpi_boundary.h \
	timestep.h \
	domain.h \
	init.h \
	refine.h \
	event.h \
	simulation.h \
	graphic.h \
	adaptive.h \
	source.h \
	tension.h \
	vof.h \
	utils.h \
	ocean.h \
	wave.h \
	levelset.h \
	isocube.h \
	cartesian.h \
	surface.h \
	unstructured.h \
	map.h \
	river.h \
	moving.h \
	balance.h \
	metric.h \
	particle.h \
	version.h

pkginclude_HEADERS = \
	$(GFS_HDS) \
	spatial.h \
	function.h

include_HEADERS = \
	gfs.h

if OPEN_MEMSTREAM
MEMSTREAM = open_memstream.c
endif

SRC = \
	ftt.c \
	fluid.c \
	variable.c \
	output.c \
	solid.c \
	poisson.c \
	advection.c \
	boundary.c \
	mpi_boundary.c \
	timestep.c \
	domain.c \
	init.c \
	refine.c \
	event.c \
	simulation.c \
	graphic.c \
	adaptive.c \
	source.c \
	tension.c \
	vof.c \
	utils.c \
	ocean.c \
	wave.c \
	levelset.c \
	myc.h \
	myc2d.h \
	cartesian.c \
	surface.c \
	unstructured.c \
	map.c \
	river.c \
	moving.c \
	balance.c \
	metric.c \
        particle.c \
	$(GFS_HDS) \
	$(MEMSTREAM)

domain.c: version.h

libgfs3D_la_LDFLAGS = $(NO_UNDEFINED)\
        -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic
libgfs3D_la_SOURCES = $(SRC)
libgfs3D_la_LIBADD = $(GTS_LIBS)

libgfs2D_la_LDFLAGS = $(NO_UNDEFINED)\
        -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic
libgfs2D_la_SOURCES = $(SRC)
libgfs2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libgfs2D_la_LIBADD = $(GTS_LIBS)

CLEANFILES = $(BUILT_SOURCES)

EXTRA_DIST = \
	ftt_internal.c \
	moving2.c \
	open_memstream.c \
	verify.h \
	m4.awk

bin_PROGRAMS = gerris2D gerris3D

gerris2D_SOURCES = gerris.c
gerris2D_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
gerris2D_LDADD = $(GFS2D_LIBS)
gerris2D_DEPENDENCIES = libgfs2D.la

gerris3D_SOURCES = gerris.c
gerris3D_LDADD = $(GFS3D_LIBS)
gerris3D_DEPENDENCIES = libgfs3D.la

if DARCS_CONTROLLED
DARCS_VERSION = darcsversion
else
DARCS_VERSION = 
endif

version.h: $(DARCS_VERSION)

darcsversion:
	sh darcsversion.sh
