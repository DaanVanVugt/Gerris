if test x$donotrun != xtrue; then
    if gerris2D $1; then :
    else
	exit 1
    fi
fi

if awk '{print $5 " & " $7 " & " $9 "\\\\"}' < norms > norms.tex && \
   awk '{print "{\\color{blue}" $5 "} & {\\color{blue}" $7 "} & {\\color{blue}" $9 "}"}' < norms.ref >> norms.tex ; then :
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('norms',3,5) - Curve('norms.ref',3,5)).max() > 0. or\
   (Curve('norms',3,7) - Curve('norms.ref',3,7)).max() > 0. or\
   (Curve('norms',3,9) - Curve('norms.ref',3,9)).max() > 0. or\
   (Curve('norms',3,11) - Curve('norms.ref',3,11)).max() > 0.:
    exit(1)
if Curve('t',1,2).max() > 2e-5:
    exit(1)
EOF
else
   exit 1
fi
