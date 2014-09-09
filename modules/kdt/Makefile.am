## Process this file with automake to produce Makefile.in

noinst_LTLIBRARIES = libkdt.la

bin_PROGRAMS = \
	xyz2kdt \
	kdt2kdt \
	kdtquery

libkdt_la_CFLAGS = $(AM_CFLAGS) -D_FILE_OFFSET_BITS=64
libkdt_la_SOURCES = \
	kdt.c \
	kdt.h

xyz2kdt_SOURCES = xyz2kdt.c kdt.h
xyz2kdt_LDADD = -lkdt -lm
xyz2kdt_CFLAGS = $(AM_CFLAGS)
xyz2kdt_DEPENDENCIES = libkdt.la

kdt2kdt_SOURCES = kdt2kdt.c kdt.h
kdt2kdt_LDADD = -lkdt -lm
kdt2kdt_CFLAGS = $(AM_CFLAGS)
kdt2kdt_DEPENDENCIES = libkdt.la

kdtquery_SOURCES = kdtquery.c kdt.h
kdtquery_LDADD = -lkdt -lm
kdtquery_CFLAGS = $(AM_CFLAGS)
kdtquery_DEPENDENCIES = libkdt.la
