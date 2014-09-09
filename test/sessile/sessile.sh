angles="30 45 60 90 120 150 180"
for level in 4 5 6; do
    for theta in $angles; do
	if gerris2D -DANGLE=$theta -DLEVEL=$level sessile.gfs; then :
	else
	    exit 1
	fi
	tail -n 1 k | awk -v theta=$theta '{
        print theta,-$5,-$7,-$11
    }'
    done > rk-$level
done

for theta in $angles; do
    rm -f convergence-$theta
    for level in 4 5 6; do
	awk -v theta=$theta -v level=$level '
        BEGIN {
          pi = 3.14159265358979
        }
        function R(theta,      vol)
        {
          vol=pi*0.3**2/2.
          return sqrt(vol/(theta-cos(theta)*sin(theta)));
        }
        {
          kexact = 1./R(theta*pi/180.)
          if ($1 == theta) print 0.3*(2**level),($3 - kexact)/kexact; 
        }' < rk-$level >> convergence-$theta
    done
done

awk 'BEGIN {
       pi = 3.14159265358979
     }
     function R(theta,      vol)
     {
       vol=pi*0.3**2/2.
       return sqrt(vol/(theta-cos(theta)*sin(theta)));
     }
     {
       theta = $1
       kexact = 1./R(theta*pi/180.)
       if ($3 > kexact) error = $3 - kexact;
       else error = kexact - $3;
       print $1,error/kexact; 
     }' < rk-6 > error-6

if gnuplot <<EOF; then :
set term postscript eps color enhanced lw 3 20
set output 'rk.eps'
set xlabel "Contact angle (degrees)"
set ylabel "Radius"
set xtics 0,30,180
vol=pi*0.3**2/2.
R(theta)=sqrt(vol/(theta-cos(theta)*sin(theta)))
plot [25:180]'rk-4' u 1:(1./\$3):(1./\$2):(1./\$4) w yerr t 'Gerris (4 levels)' pt 7, \
             'rk-5' u 1:(1./\$3):(1./\$2):(1./\$4) w yerr t 'Gerris (5 levels)' pt 9, \
             'rk-6' u 1:(1./\$3):(1./\$2):(1./\$4) w yerr t 'Gerris (6 levels)' pt 11, \
             R(x*pi/180.) t 'analytical'

set output 'convergence.eps'
set xlabel 'Spatial resolution (grid points per drop radius)'
set ylabel 'Relative error on curvature'
set logscale
set xtics 4.8,2,19.2
plot [4:25][0.0005:]\
  'convergence-30' u 1:(abs(\$2)) w lp t '30 degrees', \
  'convergence-45' u 1:(abs(\$2)) w lp t '45 degrees', \
  'convergence-60' u 1:(abs(\$2)) w lp t '60 degrees', \
  'convergence-90' u 1:(abs(\$2)) w lp t '90 degrees', \
  'convergence-120' u 1:(abs(\$2)) w lp t '120 degrees', \
  'convergence-150' u 1:(abs(\$2)) w lp lt 1 t '150 degrees', \
  'convergence-180' u 1:(abs(\$2)) w lp t '180 degrees', \
   x**-2/3. lt 0 t 'x^{-2}', x**-1/3. lt 0 t 'x^{-1}'
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('error-6',1,2) - Curve('error-6.ref',1,2)).max() > 0.:
    print (Curve('error-6',1,2) - Curve('error-6.ref',1,2)).max()
    exit(1)
EOF
else
   exit 1
fi
