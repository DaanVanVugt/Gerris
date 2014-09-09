if test x$donotrun != xtrue; then
    awk 'BEGIN{for (x = 0.35/100.; x <= 0.35; x += 0.35/100.) print x,x,0;}' > thetapi4
    cp -f electro.gfs result-7.gfs
    rm -f fprof-* convergence
    for level in 8 9 10; do
	level1=`expr $level - 1`
	echo -n "$level " >> convergence
	if sed -e 's/GfsTime.*$/Time { end = 1 }/' \
            -e "s/maxlevel = $level1/maxlevel = $level/g" < result-$level1.gfs | \
	    gerris2D - 2> log; then
	    mv -f result.gfs result-$level.gfs
	    mv -f fprof fprof-$level
	    cat log >> convergence
	else
	    cat log > /dev/stderr
	    exit 1
	fi
    done
fi

if echo "Save figure.eps { format = EPS }" | gfsview-batch2D result-10.gfs figure.gfv; then :
else
    exit 1
fi

if cat <<EOF | gnuplot; then :
set term postscript eps color enhanced lw 2 18 solid
set output 'profile.eps'

#
# Theoretical velocity profile taken from Tomar et al, 2007
# 
R=5.1
Q=10.
lambda=1.
theta=pi/4.

A=-9./10.*(R-Q)/(R+2.)**2/(1.+lambda)
vtr(x)=(x < 1. ? A*x*(1.-x**2)*(3.*sin(theta)**2-1.) : \
                 A*x**(-2)*(x**(-2)-1.)*(3.*sin(theta)**2-1.))
vtt(x)=(x < 1. ? 3.*A/2.*x*(1.-5./3.*x**2)*sin(2.*theta) : \
                 -A*x**(-4)*sin(2.*theta))

set xlabel 'r/R_{0}'
set ylabel 'v/v_{c}'
set key bottom right
set samples 1000
plot \
          'fprof-10' u 1:2 t '', \
          'fprof-10' u 1:3 t '', \
          vtr(x) t 'v_{r}', vtt(x) t 'v_{/Symbol q}'
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
c = Curve()
print (Curve('convergence',1,2) - Curve('convergence.ref',1,2)).max()
print (Curve('convergence',1,3) - Curve('convergence.ref',1,3)).max()
print Curve('rhoe',3,5).max()
if (Curve('convergence',1,2) - Curve('convergence.ref',1,2)).max() > 1e-6 or\
   (Curve('convergence',1,3) - Curve('convergence.ref',1,3)).max() > 1e-6 or\
    Curve('rhoe',3,5).max() > 1e-9:
    exit(1)
EOF
else
   exit 1
fi
