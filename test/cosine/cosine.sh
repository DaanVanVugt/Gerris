levels="4 5 6 7"
alphas="45 90"

if test x$donotrun != xtrue; then
    for alpha in $alphas; do
	for i in $levels; do
	    if gerris2D -DLEVEL=$i -DALPHA=$alpha cosine.gfs 2> log-$i-$alpha; then :
	    else
		exit 1
	    fi
	done
    done
fi

for alpha in $alphas; do
    rm -f error-$alpha
    for i in $levels; do
	tail -n 1 error-$i-$alpha >> error-$alpha
    done

    if cat <<EOF | gnuplot ; then :
set term postscript eps color enhanced lw 2 18
set output 'order-$alpha.eps'
set logscale
set xtics 16,2,256
set key spacing 1.5 bottom left
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f2(x)=a2+b2*x
fit [3:]f2(x) 'error-$alpha' u (log(2**\$1)):(log(\$4)) via a2,b2
fm(x)=am+bm*x
fit [3:]fm(x) 'error-$alpha' u (log(2**\$1)):(log(\$5)) via am,bm
set xlabel 'Spatial resolution'
set ylabel 'Relative error norms'
plot 'error-$alpha' u (2**\$1):4 t 'L2' w lp ps 2,  exp(f2(log(x))) t ftitle(a2,b2), \
     'error-$alpha' u (2**\$1):5 t 'Max' w lp ps 2, exp(fm(log(x))) t ftitle(am,bm), \
     'rossmanith$alpha' u 1:3 t 'L2 (Rossmanith)' w lp ps 2, \
     'rossmanith$alpha' u 1:4 t 'Max (Rossmanith)' w lp ps 2 lt 1

set output 'error-$alpha.eps'
unset logscale
set logscale y
set xtics auto
set key spacing 1
set xlabel 'Time'
set ylabel 'Maximum relative error'
plot 'error-4-$alpha' u 2:5 t '4 levels' w l,\
     'error-5-$alpha' u 2:5 t '5 levels' w l,\
     'error-6-$alpha' u 2:5 t '6 levels' w l,\
     'error-7-$alpha' u 2:5 t '7 levels' w l
EOF
    else
	exit 1
    fi
done

if cat <<EOF | python ; then :
from check import *
from sys import *
c = Curve()
print (Curve('error-45',1,4) - Curve('error-45.ref',1,4)).max()
print (Curve('error-45',1,5) - Curve('error-45.ref',1,5)).max()
print (Curve('error-90',1,4) - Curve('error-90.ref',1,4)).max()
print (Curve('error-90',1,5) - Curve('error-90.ref',1,5)).max() 
if (Curve('error-45',1,4) - Curve('error-45.ref',1,4)).max() > 1e-10 or\
   (Curve('error-45',1,5) - Curve('error-45.ref',1,5)).max() > 1e-10 or\
   (Curve('error-90',1,4) - Curve('error-90.ref',1,4)).max() > 1e-10 or\
   (Curve('error-90',1,5) - Curve('error-90.ref',1,5)).max() > 1e-10:
    exit(1)
EOF
else
   exit 1
fi
