if test x$donotrun != xtrue; then
    tmp=`mktemp -d`

    for level in 5 6 7; do
	if sed "s/LEVEL/$level/g" < $1 | gerris2D -; then :
            awk '{print 2**'$level', $0}' error$level.dat >> $tmp/error.dat
	else
	    exit 1
	fi
    done
    cat $tmp/error.dat > error
    rm -rf $tmp
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps enhanced color lw 3 solid 20
    set output 'divmax.eps'
    set xlabel 'Time'
    set ylabel 'Divergence Max'
    plot [0:2]'div5' u 3:9 t "5" w l, 'div6' u 3:9 t "6" w l, 'div7' u 3:9 t "7" w l, 'div5.ref' u 3:9 t "5 (ref)" w l, 'div6.ref' u 3:9 t "6 (ref)" w l, 'div7.ref' u 3:9 t "7 (ref)" w l
    set output 'divL2.eps'
    set ylabel 'Divergence L2'
    plot [0:2]'div5' u 3:7 t "5" w l, 'div6' u 3:7 t "6" w l, 'div7' u 3:7 t "7" w l, 'div5.ref' u 3:7 t "5 (ref)" w l, 'div6.ref' u 3:7 t "6 (ref)" w l, 'div7.ref' u 3:7 t "7 (ref)" w l
    set output 'kinetic.eps'
    set ylabel 'Kinetic energy'
    plot [0:2]'kinetic5' u 3:5 t "5" w l, 'kinetic6' u 3:5 t "6" w l, 'kinetic7' u 3:5 t "7" w l
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
for div in ['div5','div6','div7']:
  if (Curve(div,3,9) - Curve(div+'.ref',3,9)).max() > 0.01*Curve(div+'.ref',3,9).mean() or\
     (Curve(div,3,7) - Curve(div+'.ref',3,7)).max() > 0.01*Curve(div+'.ref',3,7).mean():
    exit(1)
EOF
else
   exit 1
fi
