## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

# change this to your OSC directory
OSC_DIR=$(HOME)/local/src/osc/home:popinet/$(PACKAGE)-snapshot

# change this to your dput command
DPUT = socksify dput

# change this to your launchpad PPA
# configured in ~/.dput.cf
PPA = gerris-ppa

# change this to the ID of the PPA GPG key
KEY = 07AB22DC

LOG = $(WEBROOT)/version/gerris_changes

SUBDIRS = src modules tools test doc desktop debian

EXTRA_DIST = ChangeLog m4 build_function_gcc build_function_libtool

if DARCS_CONTROLLED
DARCS_CHANGELOG = changelog
else
DARCS_CHANGELOG = 
endif

pkgdata_SCRIPTS = build_function libtool

ChangeLog: $(DARCS_CHANGELOG)

changelog:
	darcs changes | grep -v Ignore-this: > ChangeLog

debian-snapshot: $(DARCS_CHANGELOG) config.h src/version.h
	sh debian/changelog.sh
	dpkg-buildpackage -rfakeroot -b
	rm -f debian/repo/*
	mv -f ../*.changes ../*.deb debian/repo
	cd debian/repo && dpkg-scanpackages . /dev/null | gzip -9c > Packages.gz

launchpad-snapshot: $(DARCS_CHANGELOG) config.h src/version.h dist
	tar xzf gerris-?.?.?.tar.gz
	cd gerris-?.?.? && ./configure && sh debian/changelog.sh $(LOG) && chmod +x debian/rules && dpkg-buildpackage -rfakeroot -k$(KEY) -S -sa
	rm -f gerris-?.?.?.tar.gz
	rm -r -f gerris-?.?.?
	$(DPUT) $(PPA) gerris-snapshot_*_source.changes
	rm -f gerris-snapshot_*

buildservice-snapshot: $(DARCS_CHANGELOG) config.h src/version.h dist
	tar xzf $(PACKAGE)-?.?.?.tar.gz
	rm -f $(PACKAGE)-?.?.?.tar.gz
	mv -f $(PACKAGE)-?.?.? $(PACKAGE)-snapshot
	cd $(PACKAGE)-snapshot && ./configure && sh debian/changelog.sh $(LOG) && chmod +x debian/rules && make distclean
	tar czvf $(PACKAGE)-snapshot.tar.gz $(PACKAGE)-snapshot
	dpkg-source -b $(PACKAGE)-snapshot
	rm -r -f $(PACKAGE)-snapshot
	cd $(OSC_DIR) && osc update && rm -f $(OSC_DIR)/*
	cp -f rpm/$(PACKAGE).spec $(OSC_DIR)
	mv -f $(PACKAGE)-snapshot.tar.gz $(PACKAGE)-snapshot_*.tar.gz $(PACKAGE)-snapshot_*.dsc $(OSC_DIR)
	cd $(OSC_DIR) && osc addremove && osc commit -m "snapshot release"
