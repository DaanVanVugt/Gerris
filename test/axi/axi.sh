if test x$donotrun != xtrue; then
    rm -f error
    for level in 10 11 12 13; do
	if sed "s/LEVEL/$level/g" < $1 | \
	    gerris2D - >> error; then :
	else
	    exit 1
	fi
    done
fi

if echo "Save isolines.eps { format = EPS }" | gfsview-batch2D sim-13.gfs isolines.gfv; then :
else
    exit 1
fi

if awk '
BEGIN { n = 0 }
{
  l[n] = $1; n2[n] = $2; ni[n++] = $3;
}
END {
  for (i = 1; i < n; i++)
    print l[i] " " log(n2[i-1]/n2[i])/log(2.) " " log(ni[i-1]/ni[i])/log(2.);
}' < error > order; then :
else
    exit 1
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20
    set output 'error.eps'
    set xlabel 'Level'
    set ylabel 'Error norms'
    set key
    set logscale y
    plot 'error.ref' u 1:2 t '2 (ref)' w lp, \
         'error.ref' u 1:3 t 'max (ref)' w lp, \
         'error' u 1:2 t '2' w lp, \
         'error' u 1:3 t 'max' w lp
    set output 'order.eps'
    set xlabel 'Level'
    set ylabel 'Order'
    set key
    unset logscale
    set xtics 0,1
    set ytics 0,1
    set grid
    plot [][0:2] 'order.ref' u 1:2 t '2 (ref)' w lp, \
                 'order.ref' u 1:3 t 'max (ref)' w lp, \
                 'order' u 1:2 t '2' w lp, \
                 'order' u 1:3 t 'max' w lp
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('error',1,3) - Curve('error.ref',1,3)).max() > 1e-5:
    exit(1)
EOF
else
   exit 1
fi
