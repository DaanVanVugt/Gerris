if test x$donotrun != xtrue; then
    if gerris2D $1; then :
    else
	exit 1
    fi
fi

if python <<EOF ; then :
from check import *
from sys import *
if (Curve('xprof',3,7) - Curve('xprof.ghia',1,2)).normi() > 2e-2 or \
   (Curve('yprof',2,8) - Curve('yprof.ghia',1,2)).normi() > 1.7e-2:
    print (Curve('xprof',3,7) - Curve('xprof.ghia',1,2)).normi()
    print (Curve('yprof',2,8) - Curve('yprof.ghia',1,2)).normi()
    exit(1)
EOF
else
   exit 1
fi
