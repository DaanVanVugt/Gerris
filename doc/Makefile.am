## Process this file with automake to produce Makefile.in

SUBDIRS = tutorial examples manpages

EXTRA_DIST = figures share doxy.conf gerris.bib

noinst_PROGRAMS = autoclass

INCLUDES = -I$(top_srcdir)/src -I$(includedir) -DG_LOG_DOMAIN=\"Gfs-tools\"\
            $(GTS_CFLAGS)

autoclass: autoclass.c $(top_srcdir)/src/init.c
	$(LIBTOOL) --mode=link $(CC) $(AM_CFLAGS) $(INCLUDES) -DFTT_2D=1 \
	autoclass.c -o autoclass $(GFS2D_LIBS)

reference: Makefile.am autoclass share/autoclass.sh share/doxfilter share/renamemodules doxy.conf
	sh share/autoclass.sh | sort -k2,3 | ./autoclass > hierarchy.c
	@chmod +x share/doxfilter
	doxygen doxy.conf 2>&1 | grep -v 'explicit link request to .* could not be resolved' > doxy.log
	cd reference && sh ../share/renamemodules
	if grep [wW]arning doxy.log; then exit 1; fi

publish: reference
	tar czf reference.tgz reference
	mv -f reference.tgz $(WEBROOT)
	cd $(WEBROOT) && tar xzf reference.tgz && rm -f reference.tgz

Makefile.deps: Makefile.am ../src/Makefile.am
	ls ../src/*.[ch] | awk 'BEGIN{print "reference: \\"}{ print "\t" $$1 " \\"}' > Makefile.deps

-include Makefile.deps
