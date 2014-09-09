if test x$donotrun != xtrue; then
    for level in 2 3 4 5 6; do
	if gerris2D -DLEVEL=$level $1; then :
	else
	    exit 1
	fi
    done > error
fi

# fixme: for the moment GfsView doesn't do coordinate transformations,
# but the Gnuplot output does
if echo "Save mesh.gnu { format = Gnuplot }" | \
    gfsview-batch2D end-3.gfs mesh.gfv; then :
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
    fit f2(x) 'error' u (log(2**\$1)):(log(\$3)) via a2,b2
    fm(x)=am+bm*x
    fit fm(x) 'error' u (log(2**\$1)):(log(\$4)) via am,bm
    plot [3:96]'error' u (2**\$1):4 t 'Max' w p ps 2, exp(fm(log(x))) t ftitle(am,bm), \
               'error' u (2**\$1):3 t 'L2' w p ps 2,  exp(f2(log(x))) t ftitle(a2,b2)

    set term postscript eps lw 2 solid
    set output 'mesh.eps'
    unset key
    unset logscale
    unset border
    unset xtics
    unset ytics
    unset xlabel
    unset ylabel
    set size ratio -1
    plot 'mesh.gnu' w l lc 0
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
