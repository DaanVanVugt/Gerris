nlayers="8 16 32 64"

if test x$donotrun != xtrue; then
    rm -f error
    for nl in $nlayers; do
	if gerris2D -DNU=1. -DA=1. -DNL=$nl $1 >> error; then :
	else
	    exit 1
	fi
    done
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20 enhanced
    set output 'convergence.eps'
    set xlabel 'Spatial resolution'
    set ylabel 'Error norms'
    set logscale
    set grid
    set xtics 2
    set key spacing 1.5 bottom left
    ftitle(a,b) = sprintf("%.3g/x^{%4.2f}", exp(a), -b)
    fm(x)=am+bm*x
    fit [3:]fm(x) 'error' u (log(\$1)):(log(\$2)) via am,bm
    plot [6:80]'error' t 'Max' w p ps 2, exp(fm(log(x))) t ftitle(am,bm)
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('error',1,2) - Curve('error.ref',1,2)).max() > 1e-6:
    exit(1)
EOF
else
   exit 1
fi
