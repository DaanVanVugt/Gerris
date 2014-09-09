# catch errors
set -e

if test x$donotrun != xtrue; then
    gerris3D $1
fi

python <<EOF
from check import *
from sys import *
if (Curve('e',1,2) - Curve('e.ref',1,2)).max() > 1e-3:
    exit(1)
EOF
