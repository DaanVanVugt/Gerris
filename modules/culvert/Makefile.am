## Process this file with automake to produce Makefile.in

noinst_LTLIBRARIES = libboyd87.la

libboyd87_la_CFLAGS = $(AM_CFLAGS)
libboyd87_la_SOURCES = boyd87.c boyd87.h

noinst_PROGRAMS = check_inlet check_outlet discharge_table

check_outlet_SOURCES = check_outlet.c boyd87.h
check_outlet_LDADD = -lboyd87 -lm
check_outlet_CFLAGS = $(AM_CFLAGS)
check_outlet_DEPENDENCIES = libboyd87.la

check_inlet_SOURCES = check_inlet.c boyd87.h
check_inlet_LDADD = -lboyd87 -lm
check_inlet_CFLAGS = $(AM_CFLAGS)
check_inlet_DEPENDENCIES = libboyd87.la

discharge_table_SOURCES = discharge_table.c boyd87.h
discharge_table_LDADD = -lboyd87 -lm
discharge_table_CFLAGS = $(AM_CFLAGS)
discharge_table_DEPENDENCIES = libboyd87.la

pictures: inlet_pipe.png inlet_box.png outlet_pipe.png outlet_box.png

inlet_pipe.png: inlet_pipe.plot check_inlet
	gnuplot inlet_pipe.plot

inlet_box.png: inlet_box.plot check_inlet
	gnuplot inlet_box.plot

outlet_pipe.png: outlet_pipe.plot check_outlet
	gnuplot outlet_pipe.plot

outlet_box.png: outlet_box.plot check_outlet
	gnuplot outlet_box.plot
