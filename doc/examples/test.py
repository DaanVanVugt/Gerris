import sys
import os
import os.path
import gfs2tex

n = 0
failed = 0
for start in sys.argv[1:]:
    for root, dirs, files in os.walk(start,topdown=True):
        if not ".xvpics" in root:
            example = gfs2tex.Example(root)
            status,msg = example.test()
            if status != None:
                print "FAIL:",root
                if len(msg) > 0:
                    print " ".join(msg)
                failed += 1
            else:
                print "PASS:",root
            n += 1

if failed:
    msg = repr(failed) + " of " + repr(n) + " tests failed"
else:
    msg = "All " + repr(n) + " tests passed"

print len(msg)*"="
print msg
print len(msg)*"="

if failed:
    sys.exit(1)
