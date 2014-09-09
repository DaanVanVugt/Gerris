## Process this file with automake to produce Makefile.in

EXTRA_DIST = tutorial.tex depend.awk

Makefile.deps: Makefile.am tutorial.tex
	awk -f depend.awk file="tutorial" < tutorial.tex > Makefile.deps

clean-generic:
	$(RM) *.dvi *.aux *.log *.toc *.out \
	boundaries.pdf direction.eps direction.pdf dxscreen.eps dxscreen.pdf \
	gfs2oogl.pdf half-cylinder.eps half-cylinder.pdf refined1.pdf refined1_cells.pdf \
	refined2.pdf vorticity.eps vorticity.pdf gfsview.eps tutorial1.tex \
	Makefile.deps

tutorial.tar.gz: tutorial1.dvi
	rm -r -f tutorial
	sed 's/input{pdf.tex}/usepackage{graphicx}\\newcommand{\\gfx}{eps}/g' < tutorial.tex | sed "s/GFS_VERSION/`$(top_srcdir)/src/gerris2D -V 2>&1 | awk '{ if ($$5 == "version") print $$6}'`/g" | sed 's/\\today/'"`date +\"%B %e, %Y\"`/g" > tutorial1.tex
	hevea -fix tutorial1.tex
	imagen -res 600 -extra "pnmscale 0.24" tutorial1
	mkdir tutorial
	mv -f tutorial1[0-9][0-9][0-9].png tutorial
##	fixme: the character conversion below is a workaround for a bug in hevea version < 1.09
	konwert iso1-utf8 < tutorial1.html > tutorial/tutorial1.html
	cp ../share/darcs.css tutorial/tutorial.css
	rm -f tutorial1.h{tml,aux,ind,toc} tutorial1.image.tex
	tar cf tutorial.tar tutorial
	gzip -f --best tutorial.tar

tutorial1.dvi: tutorial.tex Makefile.deps
	sed "s/GFS_VERSION/`$(top_srcdir)/src/gerris2D -V 2>&1 | awk '{ if ($$5 == "version") print $$6}'`/g" < tutorial.tex | sed 's/\\today/'"`date +\"%B %e, %Y\"`/g" > tutorial1.tex
	latex -interaction=nonstopmode tutorial1.tex > /dev/null 2>&1
	latex -interaction=nonstopmode tutorial1.tex > /dev/null 2>&1
	latex -interaction=nonstopmode tutorial1.tex

tutorial.pdf: tutorial1.dvi
	dvips -Ppdf -G0 tutorial1.dvi -o tutorial1.ps
	ps2pdf -sPAPERSIZE=a4 -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true tutorial1.ps tutorial.pdf
	rm -f tutorial1.ps

.xmgr.eps:
	/usr/local/grace/bin/xmgrace -noask -hardcopy -hdevice EPS -printfile $@ $<

.agr.eps:
	/usr/local/grace/bin/xmgrace -noask -hardcopy -hdevice EPS -printfile $@ $<

.fig.eps:
	fig2dev -L ps $< > $@ ; sh epsbbox.sh $@ 

.jpeg.eps:
	convert $< $@.gif
	convert $@.gif EPS:$@
	- rm -f $@.gif

.jpg.eps:
	convert $< $@.gif
	convert $@.gif EPS:$@
	- rm -f $@.gif

.gif.eps:
	convert $< EPS:$@

.png.eps:
	convert $< EPS:$@

.jpg.pdf:
	convert $< EPDF:$@

.png.pdf:
	convert $< EPDF:$@

.epsf.pdf:
	epstopdf $<

.eps.pdf:
	epstopdf $<

-include Makefile.deps
