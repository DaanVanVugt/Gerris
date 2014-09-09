## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = parallel-tests

EXTRA_DIST = \
	template.tex \
	depend.py \
	test.py \
	check.py \
	gfs2tex \
	summary.sh \
	Makefile.deps

clean-generic:
	mv -f summary.sh summary.sh.bak
	$(RM) *.dvi *.aux *.log *.toc *.out tests.tex *.pyc *.bbl *.blg $(TESTS) Makefile.deps
	mv -f summary.sh.bak summary.sh

DOC = tests

# we need this for gfs.sty
TEXINPUTS = .:$(top_srcdir)/doc/examples:
# and this for bibtex
BIBINPUTS = $(top_srcdir)/doc

Makefile.deps: template.tex depend.py Makefile.am
	python depend.py `grep \\test{[^/]*} < template.tex | sed 's/\\\\test{\(.*\)}/\1/g'` > Makefile.deps
	echo "" >> Makefile.deps
	echo -n "TESTDIRS = " >> Makefile.deps
	echo `grep \\test{[^/]*} < template.tex | sed 's/\\\\test{\(.*\)}/\1/g'` >> Makefile.deps

include Makefile.deps

tests.tex: template.tex Makefile.deps
	rm -r -f tests
	python gfs2tex $(TESTDIRS)
	sed -e 's/\\test{\(.*\)\/\(.*\)}/\\input{\1\/\2\/\2.tex}/g' \
	    -e 's/\\test{\(.*\)}/\\input{\1\/\1.tex}/g' < template.tex > tests.tex

tests: tests.dvi
	if ! test -d $(DOC); then mkdir $(DOC); fi
	TEXINPUTS=$(TEXINPUTS) hevea -fix -I $(top_srcdir)/doc/examples $(DOC).tex
	TEXINPUTS=$(TEXINPUTS) imagen -res 500 -extra "pnmscale 0.24" $(DOC)
	hacha $(DOC).html
	rm -f $(DOC).html
	mv -f $(DOC)[0-9][0-9][0-9].png $(DOC)
##	fixme: the character conversion below is a workaround for a bug in hevea version < 1.09
	for f in *.html; do konwert iso1-utf8 < $$f > $(DOC)/$$f; rm -f $$f; done
	cat $(DOC).css ../doc/share/darcs.css > $(DOC)/$(DOC).css
	sh ../doc/share/fixnav.sh $(DOC)
	cp -f ../doc/share/contents.png ../doc/share/next.png ../doc/share/prev.png $(DOC)
	rm -f *_motif.gif $(DOC).h{tml,aux,ind,toc} $(DOC).image.tex $(DOC).css
	sh ../doc/examples/crossref.sh --url=http://gerris.dalembert.upmc.fr/gerris/tests/tests $(TESTDIRS)
	mv references tests

tests.dvi: tests.tex $(top_srcdir)/doc/gerris.bib
	- TEXINPUTS=$(TEXINPUTS) latex -interaction=nonstopmode tests.tex > /dev/null 2>&1
	BIBINPUTS=$(BIBINPUTS) bibtex tests
	- TEXINPUTS=$(TEXINPUTS) latex -interaction=nonstopmode tests.tex > /dev/null 2>&1
	TEXINPUTS=$(TEXINPUTS) latex -interaction=nonstopmode tests.tex

tests.pdf: tests.dvi
	dvips -Ppdf -G0 tests.dvi -o tests.ps
	ps2pdf -sPAPERSIZE=a4 -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true tests.ps tests.pdf
	rm -f tests.ps

tests.tar.gz: tests $(DOCS)
	tar czf tests.tar.gz tests $(DOCS)
	rm -r -f tests

publish: tests.tar.gz
	cp tests.tar.gz $(WEBROOT)/tests
	cd $(WEBROOT)/tests && tar xzf tests.tar.gz && rm -f tests.tar.gz
