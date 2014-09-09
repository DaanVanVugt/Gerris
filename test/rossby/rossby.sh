levels="4 5"

if test x$donotrun != xtrue; then
    cp -f eh-6.ref eh-6
    for level in $levels; do
	rm -f eh-$level
	if gerris2D -DLEVEL=$level $1; then :
	else
	    exit 1
	fi
    done
fi

for level in $levels; do
    if gnuplot <<EOF; then :
set term postscript eps lw 2 18 color

set output 'ehpm-$level.eps'
unset key
set size ratio -1
set xtics -90,45,90
set ytics -90,45,90
plot [-90:90][-90:90]'ehp-$level.gnu' w l, 'ehm-$level.gnu' w l

set output 'h-$level.eps'
plot [-90:90][-90:90]'h-$level.gnu' w l, 'href-$level.gnu' w l
EOF
    else
	exit 1
    fi
done

for i in 4 5 6; do
    echo -n $i" "
    tail -n 1 eh-$i | awk '{ print $3,$5,$7,$9 }'
done > error

if gnuplot <<EOF; then :
set term postscript eps lw 2 18 color enhanced

set output 'ec.eps'
set xlabel 'Time (days)'
set ylabel 'Total kinetic energy'
set key bottom left
plot [0:24][0:] 'ec-6' u (\$3/86400.):5 w l t '6 levels', 'ec-5' u (\$3/86400.):5 w l t '5 levels', 'ec-4' u (\$3/86400.):5 w l t '4 levels'

set output 'eh.eps'
set key top right
set ylabel 'Maximum relative error on height'
set logscale y
set grid
plot [0:24][1e-3:] 'eh-6' u (\$3/86400.):9 w l t '6 levels', 'eh-5' u (\$3/86400.):9 w l t '5 levels', 'eh-4' u (\$3/86400.):9 w l t '4 levels'

set output 'order.eps'
set logscale
set xtics 8,2,128
set key spacing 1.5 bottom left
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f2(x)=a2+b2*x
fit f2(x) 'error' u (log(2**\$1)):(log(\$4)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'error' u (log(2**\$1)):(log(\$5)) via am,bm
set xlabel 'Spatial resolution'
set ylabel 'Relative error norms'
plot 'error' u (2**\$1):5 t 'Max' w lp ps 2, exp(fm(log(x))) t ftitle(am,bm), \
     'error' u (2**\$1):4 t 'L2' w lp ps 2,  exp(f2(log(x))) t ftitle(a2,b2)
EOF
else
    exit 1
fi

if python <<EOF ; then :
from check import *
from sys import *
if (Curve('eh-4',3,9) - Curve('eh-4.ref',3,9)).max() > 1e-5 or\
   (Curve('eh-5',3,9) - Curve('eh-5.ref',3,9)).max() > 1e-5:
    print (Curve('eh-4',3,9) - Curve('eh-4.ref',3,9)).max()
    print (Curve('eh-5',3,9) - Curve('eh-5.ref',3,9)).max()
    exit(1)
EOF
else
   exit 1
fi
