#!/bin/sh

if awk 'BEGIN {
  for (x = -4.1; x <= 4.; x += 8./128.)
    for (y = -4.1; y <= 4.; y += 8./128.) {
      r = sqrt (x*x + y*y);
      print x,y,r*r/8.+cos(3.14159265356*r)/2.;
    }      
}' | xyz2kdt terrain; then :
else
    exit 1
fi

if awk 'BEGIN {
  for (x = -4.1; x <= 4.; x += 8./512.)
    for (y = -4.1; y <= 4.; y += 8./512.) {
      r = sqrt (x*x + y*y);
      if (r < 2.)
        print x,y,r*r/8.+cos(3.14159265356*r)/2.;
    }      
}' | xyz2kdt terrain-high; then :
else
    exit 1
fi

rm -f error-t error-h
for level in 4 5 6 7; do 
    if gerris2D -DLEVEL=$level terrain.gfs; then :
    else
	exit 1
    fi
done

if gnuplot <<EOF; then :
set term postscript eps color enhanced lw 2 18
set output 'error-t.eps'
set logscale
set xtics 8,2,256
set key spacing 1.5 bottom left
ftitle(a,b) = sprintf("x^{%4.2f}", b)
f2(x)=a2+b2*x
fit [:4.5]f2(x) 'error-t' u (log(2**\$1)):(log(\$3)) via a2,b2
fm(x)=am+bm*x
fit [:4.5]fm(x) 'error-t' u (log(2**\$1)):(log(\$4)) via am,bm
set xlabel 'Spatial resolution'
set ylabel 'Relative error norms'
plot [][1e-3:]'error-t' u (2**\$1):3 t 'L2' w lp ps 2,  exp(f2(log(x))) t ftitle(a2,b2), \
     'error-t' u (2**\$1):4 t 'Max' w lp ps 2, exp(fm(log(x))) t ftitle(am,bm)

set output 'error-h.eps'
fit [:4.5]f2(x) 'error-h' u (log(2**\$1)):(log(\$3)) via a2,b2
fm(x)=am+bm*x
fit [:4.5]fm(x) 'error-h' u (log(2**\$1)):(log(\$4)) via am,bm
set xlabel 'Spatial resolution'
set ylabel 'Relative error norms'
plot [][1e-3:]'error-h' u (2**\$1):3 t 'L2' w lp ps 2,  exp(f2(log(x))) t ftitle(a2,b2), \
     'error-h' u (2**\$1):4 t 'Max' w lp ps 2, exp(fm(log(x))) t ftitle(am,bm)
EOF
else
    exit 1
fi

if python <<EOF ; then :
from check import *
from sys import *
c = Curve()
if (Curve('error-t',1,3) - Curve('error-t.ref',1,3)).max() > 0. or\
   (Curve('error-t',1,4) - Curve('error-t.ref',1,4)).max() > 0. or\
   (Curve('error-h',1,3) - Curve('error-h.ref',1,3)).max() > 0. or\
   (Curve('error-h',1,4) - Curve('error-h.ref',1,4)).max() > 0.:
    print (Curve('error-t',1,3) - Curve('error-t.ref',1,3)).max()
    print (Curve('error-t',1,4) - Curve('error-t.ref',1,4)).max()
    print (Curve('error-h',1,3) - Curve('error-h.ref',1,3)).max()
    print (Curve('error-h',1,4) - Curve('error-h.ref',1,4)).max() 
    exit(1)
EOF
else
   exit 1
fi
