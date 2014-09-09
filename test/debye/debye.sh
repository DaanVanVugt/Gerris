if test x$donotrun != xtrue; then
    if gerris2D debye.gfs; then :
    else
	exit 1
    fi
fi

if cat <<EOF | gnuplot; then :
set term postscript eps color enhanced lw 2 18 solid
set output 'profile.eps'
 
set xlabel 'x'
set ylabel ' '
set key top right
plot   'profile' u 1:2 t '', 'analytical' u 1:2 t '{/Symbol f}' w lines, \
       'profile' u 1:3 t '', 'analytical' u 1:3 t '{n_{+}}' w lines, \
       'profile' u 1:4 t '', 'analytical' u 1:4 t '{n_{-}}' w lines lt 1

EOF
else
    exit 1
fi
 
if python <<EOF ; then :
from check import *
from sys import *
print (Curve('profile',1,2) - Curve('analytical',1,2)).normi()
print (Curve('profile',1,3) - Curve('analytical',1,3)).normi()
print (Curve('profile',1,4) - Curve('analytical',1,4)).normi()
if (Curve('profile',1,2) - Curve('analytical',1,2)).normi() > 1.6e-3 or \
   (Curve('profile',1,3) - Curve('analytical',1,3)).normi() > 6.9e-3 or \
   (Curve('profile',1,4) - Curve('analytical',1,4)).normi() > 6.77e-3:
    exit(1)
EOF
else
   exit 1
fi
