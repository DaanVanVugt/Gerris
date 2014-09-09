emax=$2

if test x$donotrun != xtrue; then
    awk '{printf("%4.3f %4.3f 0.\n",0.,$1)}' < prof.ref > profile
    if ( gerris2D -DLEVEL=7 $1 ) ; then :
    else
	exit 1
    fi     	
    
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
avg = (Curve('prof.ref',1,2) - Curve('prof.dat',3,5)).mean()
print avg
c = Curve()
for p in Curve('prof.dat',3,5).l:
    c.l.append((p[0], p[1] + avg))
if (Curve('prof.ref',1,2) - c).normi() > $emax:
   print (Curve('prof.ref',1,2) - c).normi()
   exit(1)
EOF
else
   exit 1
fi
