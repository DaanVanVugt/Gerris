## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/src -I$(includedir) -DG_LOG_DOMAIN=\"Gfs-tools\"\
            $(GTS_CFLAGS)

bin_PROGRAMS = \
	gfs2oogl2D \
	gfs2oogl3D \
	streamanime \
	ppmcombine \
	gfscompare2D \
	gfscompare3D \
	gfscombine2D \
	gfscombine3D \
	gfsjoin2D \
	gfsjoin3D \
	shapes

bin_SCRIPTS = \
	darcs2dist \
	bat2gts \
	ppm2mpeg \
	ppm2theora \
	ppm2video \
	gfs2gfs \
	gfsjoin \
	gfsplot

EXTRA_DIST = $(bin_SCRIPTS)

gfs2oogl2D_SOURCES = gfs2oogl.c
gfs2oogl3D_SOURCES = gfs2oogl.c
gfscompare2D_SOURCES = gfscompare.c
gfscompare3D_SOURCES = gfscompare.c
gfscombine2D_SOURCES = gfscombine.c
gfscombine3D_SOURCES = gfscombine.c
gfsjoin2D_SOURCES = gfsjoin2.c
gfsjoin3D_SOURCES = gfsjoin2.c

gfs2oogl2D_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
gfs2oogl2D_LDADD = $(GFS2D_LIBS)
gfs2oogl3D_LDADD = $(GFS3D_LIBS)
gfscompare2D_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
gfscompare2D_LDADD = $(GFS2D_LIBS)
gfscompare3D_LDADD = $(GFS3D_LIBS)
gfscombine2D_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
gfscombine2D_LDADD = $(GFS2D_LIBS)
gfscombine3D_LDADD = $(GFS3D_LIBS)
gfsjoin2D_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
gfsjoin2D_LDADD = $(GFS2D_LIBS)
gfsjoin3D_LDADD = $(GFS3D_LIBS)

ppmcombine_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
ppmcombine_LDADD = $(GFS2D_LIBS)

streamanime_LDADD = $(GFS3D_LIBS)

shapes_LDADD = $(GTS_LIBS)
