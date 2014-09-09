#!/bin/sh

rm -f error
for level in 5 6 7; do
    if gerris2D -DLEVEL=$level $1 >> error; then :
    else
	exit 1
    fi
done

if cat <<EOF | gnuplot ; then :
    set term postscript eps enhanced color lw 3 solid 20
    set output 'kinetic.eps'
    set xlabel 'Time'
    set ylabel 'Kinetic energy'
    plot 'kinetic-5' u 3:5 t "5" w l, 'kinetic-6' u 3:5 t "6" w l, 'kinetic-7' u 3:5 t "7" w l
    set output 'accuracy.eps'
    set logscale 
    set ylabel 'Relative error norms'
    set xlabel 'Spatial resolution'
    ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
    f2(x)=a2+b2*x
    fit f2(x) 'error' u (log(\$1)):(log(\$4)) via a2,b2
    fm(x)=am+bm*x
    fit fm(x) 'error' u (log(\$1)):(log(\$5)) via am,bm
    set xrange[25:150]
    set xtics 32,2,128
    set key spacing 1.5 top right
    plot 'error' u (\$1):4 t 'L2' w p ps 2, exp(f2(log(x))) t ftitle(a2,b2), \
         'error' u (\$1):5 t 'Lmax' w p ps 2, exp(fm(log(x))) t ftitle(am,bm)
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('error',1,4) - Curve('error.ref',1,4)).max() > 0.01*Curve('error.ref',1,4).mean() or\
   (Curve('error',1,5) - Curve('error.ref',1,5)).max() > 0.01*Curve('error.ref',1,5).mean():
  print (Curve('error',1,4) - Curve('error.ref',1,4)).max(), 0.01*Curve('error.ref',1,4).mean()
  print (Curve('error',1,5) - Curve('error.ref',1,5)).max(), 0.01*Curve('error.ref',1,5).mean()
  exit(1)
EOF
else
   exit 1
fi

