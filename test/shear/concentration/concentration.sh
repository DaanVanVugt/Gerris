for level in 5 6 7; do
    if gerris2D -DLEVEL=$level concentration.gfs; then :
    else
	exit 1
    fi
done

for i in 1 2; do
    echo "Save t$i-half.eps { format = EPS }" | gfsview-batch2D half-7.gfs t$i.gfv
    echo "Save t$i-end.eps { format = EPS }" | gfsview-batch2D end-7.gfs t$i.gfv
done

rm -f convergence convergence1 convergence2
for level in 5 6 7; do
    for i in "" 1 2; do
	awk -v level=$level '{ print level, $5, $7, $9}' < dt$i-$level >> convergence$i
    done
done

if gnuplot <<EOF; then :
set term postscript eps color lw 2 18 solid
set output 'convergence.eps'
set logscale y
set xlabel 'Level'
set ylabel 'Error'
set xtics 4,1,8
set key bottom left reverse Left
set pointsize 1.5
plot [4.5:7.5][1e-5:]\
     'convergence'  u 1:2 t 'T (L1)', 'convergence'  u 1:4 t 'T (max)', \
     'convergence1' u 1:2 t 'C (L1)', 'convergence1' u 1:4 t 'C (max)', \
     'convergence2' u 1:2 t 'G (L1)', 'convergence2' u 1:4 t 'G (max)' lc 1, \
     30./2**x t 'first-order', \
     6./4**x t 'second-order'
EOF
else
    exit 1
fi

for i in "" 1 2; do
    if python <<EOF ; then :
from check import *
from sys import *
if (Curve('convergence$i',1,2) - Curve('convergence$i.ref',1,2)).max() > 0. or\
   (Curve('convergence$i',1,3) - Curve('convergence$i.ref',1,3)).max() > 0. or\
   (Curve('convergence$i',1,4) - Curve('convergence$i.ref',1,4)).max() > 0.:
  exit(1)
EOF
    else
	exit 1
    fi
done
