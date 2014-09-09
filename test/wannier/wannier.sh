if test x$donotrun != xtrue; then
    for level in 4 5 6 7 8; do
	if gerris2D -DLEVEL=$level $1; then :
	else
	    exit 1
	fi
    done > error
fi

if echo "Save streamlines.eps { format = EPS line_width = 0.5 }" | \
    gfsview-batch2D end-8.gfs streamlines.gfv; then :
else
    exit 1
fi

if gnuplot <<EOF ; then :
    set term postscript eps color lw 3 solid 20 enhanced
    set output 'convergence.eps'
    set xlabel 'Spatial resolution'
    set ylabel 'Error norms'
    set logscale
    set grid
    set xtics 2
    set key spacing 1.5 bottom left
    ftitle(a,b) = sprintf("%.3g/x^{%4.2f}", exp(a), -b)
    f2(x)=a2+b2*x
    fit [3:]f2(x) 'error' u (log(2**\$1)):(log(\$3)) via a2,b2
    fm(x)=am+bm*x
    fit [3:]fm(x) 'error' u (log(2**\$1)):(log(\$4)) via am,bm
    plot [10:384]'error' u (2**\$1):4 t 'Max' w p ps 2, exp(fm(log(x))) t ftitle(am,bm), \
                 'error' u (2**\$1):3 t 'L2' w p ps 2,  exp(f2(log(x))) t ftitle(a2,b2)
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('error',1,4) - Curve('error.ref',1,4)).max() > 1e-6:
    exit(1)
EOF
else
   exit 1
fi
