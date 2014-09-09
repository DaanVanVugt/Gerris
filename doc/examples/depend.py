#!/usr/bin/python

import sys
import os
import os.path
import gfs2tex

print "## File generated automatically by depend.py: do not modify by hand"
print ""

dists = ""
depends = ""
docs = ""
tests = ""
for start in sys.argv[1:]:
    tests += "\\\n\t" + start + ".sh"
    f = open(start + ".sh", "w")
    f.write("python -u test.py " + start + "\n")
    f.close()
    os.chmod(start + ".sh",0755)
    for root, dirs, files in os.walk(start,topdown=True):
        if not ".xvpics" in root:
            example = gfs2tex.Example(root)
            name = example.path + "/" + example.name + ".gfs"
            docs += "\\\n\t" + name + ".html"
            dists += "\\\n\t" + name
            depends += "\\\n\t" + name
            for f in example.required:
                if f != "S60-scaled.gts":
                    dists += "\\\n\t" + example.path + "/" + f
            for f in example.generated:
                depends += "\\\n\t" + example.path + "/" + f
                if f[-4:] == ".mpg" or f[-4:] == ".ogv" or f[-4:] == ".png" or f[-4:] == ".mp4":
                    docs += "\\\n\t" + example.path + "/" + f

print "DOCS = " + docs + dists
print ""
print "TESTS = " + tests
print ""
print "EXTRA_DIST += $(TESTS)" + dists
print ""
print "examples.tex: " + depends
