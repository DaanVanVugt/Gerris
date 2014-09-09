if test x$donotrun != xtrue; then
    if gerris2D $1; then :
    else
	exit 1
    fi
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
print (Curve('xprof',3,7) - Curve('../xprof.ghia',1,2)).normi()
print (Curve('yprof',2,8) - Curve('../yprof.ghia',1,2)).normi()
if (Curve('xprof',3,7) - Curve('../xprof.ghia',1,2)).normi() > 2.2e-2 or \
   (Curve('yprof',2,8) - Curve('../yprof.ghia',1,2)).normi() > 2.1e-2:
    exit(1)
EOF
else
   exit 1
fi
