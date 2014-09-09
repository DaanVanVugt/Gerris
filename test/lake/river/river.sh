ratio=1
if test x$donotrun != xtrue; then
    for nl in 4 8 16 32; do
	if gerris2D -DRATIO=$ratio -DNL=$nl -DRE=1 -DLEVEL=5 river.gfs; then :
	else
	    exit 1
	fi
    done
fi

for nl in  4 8 16 32; do
    awk -v nl=$nl '{
      z = $3
      u = z/4.*(3.*z-2.)
      e = $4 - u
      if (e < 0.) e = -e;
      if (e > emax) emax = e;
    } END {
      print nl,emax
    }' < uprof-$ratio-1-$nl
done > error

if gnuplot <<EOF; then :
set term postscript eps color lw 2 18 solid enhanced
set output 'uprof.eps'
set key top left
set xlabel 'z'
set ylabel 'u'
plot [0:1]x/4.*(3.*x-2.) t 'analytical', \
          'uprof-$ratio-1-4' u 3:4 t '4 layers', \
          'uprof-$ratio-1-8' u 3:4 t '8 layers', \
          'uprof-$ratio-1-16' u 3:4 t '16 layers', \
          'uprof-$ratio-1-32' u 3:4 t '32 layers'

set output 'error.eps'
fm(x)=am+bm*x
fit fm(x) 'error' u (log(\$1)):(log(\$2)) via am,bm
set key spacing 1.5 top right
ftitle(a,b) = sprintf("%.2g/x^{%4.2f}", exp(a), -b)
set xlabel 'Spatial resolution'
set ylabel 'Maximum error'
set logscale
set xtics 2,2,64
plot [2:64]'error' u 1:2 t '' w lp ps 2, exp(fm(log(x))) t ftitle(am,bm)
EOF
else
    exit 1
fi

if python <<EOF ; then :
from check import *
from sys import *
if (Curve('error',1,2) - Curve('error.ref',1,2)).max() > 1e-6:
    print (Curve('error',1,2) - Curve('error.ref',1,2)).max()
    exit(1)
EOF
else
   exit 1
fi
