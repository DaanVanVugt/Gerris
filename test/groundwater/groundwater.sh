if test x$donotrun != xtrue; then
    rm -f p U V
    for L in `seq 9`; do
	if gerris2D -DLEVEL=$L $1 ; then :
	else
	    exit 1
	fi
    done
    if echo Save solution.eps '{ format = EPS }' | \
	gfsview-batch2D solution.gfs `basename $1 .gfs`.gfv ; then :
    else
	exit 1
    fi
else
    exit 1
fi

if cat <<EOF | gnuplot ; then :
    set xlabel 'Refinement level'
    set ylabel 'Error norm'
    set logscale y
    set terminal postscript eps color lw 2 solid 18 enhanced
    set output 'convergence.eps'
    set style data lp
    set key bottom left
    plot 'p' t '1 (p)', \
        '' u 1:4 t 'max (p)', \
        'U' t '1 (U)' w l, \
        '' u 1:4 t 'max (U)' w l, \
        'V' t '1 (V)' w p, \
        '' u 1:4 t 'max (V)' w p lt 2, \
        5./2**x t '1/x', \
        5e-2/2**(2*x) t '1/x^2'
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('p',1,2) - Curve('p.ref',1,2)).max() > 1e-8 or\
   (Curve('p',1,3) - Curve('p.ref',1,3)).max() > 1e-6:
    print (Curve('p',1,2) - Curve('p.ref',1,2)).max()
    print (Curve('p',1,3) - Curve('p.ref',1,3)).max()
    exit(1)
if (Curve('U',1,2) - Curve('U.ref',1,2)).max() > 1e-8 or\
   (Curve('U',1,3) - Curve('U.ref',1,3)).max() > 1e-6:
    print (Curve('U',1,2) - Curve('U.ref',1,2)).max()
    exit(1)
EOF
else
    exit 1
fi
