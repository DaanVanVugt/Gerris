levels="5 6 7 8"

rm -f error laplace

# try several small perturbations of the interface position
# the results are quite sensitive to interface configurations
# for i in 4 3 2 1 0; do
for i in 0; do
    diameter=`awk -v i=$i 'BEGIN { print 0.2 - i/512.}'`
    rm -f fit
    if test x$donotrun != xtrue; then
	for level in $levels; do
	    if gerris2D -D LEVEL=$level -D DIAMETER=$diameter oscillation.gfs >> fit; then :
	    else
		exit 1
	    fi
	done
    fi

    if awk -v D=$diameter 'BEGIN {
              n = 2.
              sigma = 1.
              rhol = 1.
              rhog = 1./1000.
              r0 = D/2.
              omega0 = sqrt((n**3-n)*sigma/((rhol+rhog)*r0**3))
              empirical_constant = 30.
            }{
              print D*2.**$1, $4/2./omega0-1., D >> "error"
              print D*2.**$1, (1./($3**2.*D**3.))*empirical_constant**2, D >> "laplace"
            }' < fit; then :
    else
	exit 1
    fi
done

rm -f fit-*
if awk '{
    level = $1; a = $2; b = $3; c = $4;
    for (t = 0; t <= 1.; t += 0.005)
      print t, 2.*a*exp(-b*t) >> "fit-" level;
}' < fit ; then :
else
    exit 1
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20 enhanced

    set output 'k.eps'
    set xlabel 'Time'
    set ylabel 'Kinetic energy'
    set logscale y
    plot [0:1][8e-5:]'k-8' u 3:5 t "256x256" w l, 'k-7' u 3:5 t "128x128" w l, 'k-6' u 3:5 t "64x64" w l, 'k-5' u 3:5 t "32x32" w l, 'fit-8' t "fit" w l lt 7, 'fit-7' t "" w l lt 7, 'fit-6' t "" w l lt 7, 'fit-5' t "" w l lt 7

    set output 'laplace.eps'
    set xlabel 'Diameter (grid points)'
    set ylabel 'Equivalent Laplace number'
    set logscale y
    set logscale x 2
    set grid
    plot 'laplace' t "" w p pt 5 ps 2

    set output 'frequency.eps'
    set xlabel 'Diameter (grid points)'
    set ylabel 'Frequency error (%)'
    unset grid
    set xzeroaxis
    set key spacing 1.5 top right
    ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
    f(x)=a+b*x
    fit f(x) 'error' u (log(\$1)):(log(abs(\$2)*100.)) via a,b
    plot 'error' u (\$1):(abs(\$2)*100.) t "" w p pt 5 ps 2, \
         exp(f(log(x))) t ftitle(a,b)
    
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('fit',1,3) - Curve('fit.ref',1,3)).max() > 1e-2 or\
   (Curve('fit',1,4) - Curve('fit.ref',1,4)).max() > 1e-2:
  exit(1)
EOF
else
   exit 1
fi
