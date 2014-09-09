#!/bin/bash

# exit on errors
set -e

if test x$donotrun != xtrue; then
    for omega in 0 1; do
	gerris2D -DOMEGA=$omega nonlinear.gfs
	gerris2D -DOMEGA=$omega river.gfs
	gerris3D -DOMEGA=$omega ocean.gfs
    done
fi

for  omega in 0 1; do
    gnuplot <<EOF
set term postscript eps color solid lw 2 18
set output 'error-$omega.eps'
set xlabel 'Time'
set ylabel 'Maximum relative error'
set xrange [0.01:]
set logscale y
plot 'error-$omega' u 3:9 w l t 'Incompressible', \
     'error-ocean-$omega' u 3:9 w l t 'Linearised free surface', \
     'error-river-$omega' u 3:9 w l t 'Saint-Venant'
EOF

    echo "Save end-$omega.eps { format = EPS }" | gfsview-batch2D end-$omega.gfs error.gfv
    echo "Save end-river-$omega.eps { format = EPS }" | \
	gfsview-batch2D end-river-$omega.gfs error.gfv
    echo "Save end-ocean-$omega.eps { format = EPS }" | \
	gfsview-batch3D end-ocean-$omega.gfs error-ocean.gfv
done

python <<EOF
from check import *
from sys import *
if (Curve('error-1',3,9) - Curve('error-1.ref',3,9)).max() > 1e-7 or\
   (Curve('error-river-1',3,9) - Curve('error-river-1.ref',3,9)).max() > 1e-7 or\
   (Curve('error-ocean-1',3,9) - Curve('error-ocean-1.ref',3,9)).max() > 1e-7:
    print (Curve('error-1',3,9) - Curve('error-1.ref',3,9)).max()
    print (Curve('error-river-1',3,9) - Curve('error-river-1.ref',3,9)).max()
    print (Curve('error-ocean-1',3,9) - Curve('error-ocean-1.ref',3,9)).max()
    exit(1)
EOF
