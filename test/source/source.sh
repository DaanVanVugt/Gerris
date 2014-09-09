levels="7 8 9 10 11"

if test x$donotrun != xtrue; then
    for i in $levels; do
	if gerris2D -DLEVEL=$i source.gfs 2> log-$i; then :
	else
	    exit 1
	fi
    done
fi

rm -rf error
for i in $levels; do
    tail -n 1 error-$i >> error
done

if gnuplot <<EOF ; then :
set term postscript eps color enhanced lw 3 solid 20
set output 'error.eps'
set logscale
set xlabel 'Spatial resolution'
set ylabel 'Relative error norms'
set xtics 128,2,2048
set key spacing 1.5 bottom left
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f2(x)=a2+b2*x
fit [5:]f2(x) 'error' u (log(2**\$1)):(log(\$3)) via a2,b2
fm(x)=am+bm*x
fit [5:]fm(x) 'error' u (log(2**\$1)):(log(\$4)) via am,bm
plot 'error' u (2**\$1):3 t 'L2' w p ps 2, exp(f2(log(x))) t ftitle(a2,b2), \
     'error' u (2**\$1):4 t 'Lmax' w p ps 2, exp(fm(log(x))) t ftitle(am,bm)
EOF
else
: #    exit 1
fi

if echo "Save stdout { width = 800 height = 800 }" | \
    gfsview-batch2D end-11.gfs source.gfv | \
    convert -colors 256 ppm:- velfield.eps; then :
else
    exit 1
fi

if echo "Save stdout { width = 800 height = 800 }" | \
    gfsview-batch2D end-11.gfs error.gfv | \
    convert -colors 256 ppm:- localerror.eps; then :
else
    exit 1
fi

if python <<EOF ; then :
from check import *
from sys import *
if (Curve('error',1,4) - Curve('error.ref',1,4)).max() > 1e-4 or\
   (Curve('error',1,3) - Curve('error.ref',1,3)).max() > 1e-4:
    print (Curve('error',1,4) - Curve('error.ref',1,4)).max()
    print (Curve('error',1,3) - Curve('error.ref',1,3)).max()
    exit(1)
EOF
else
   exit 1
fi
