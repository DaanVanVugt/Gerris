if test x$donotrun != xtrue; then
    for level in 5 6 7 8; do
	if gerris2D -DLEVEL=$level $1; then :
	else
	    exit 1
	fi
    done
fi

for level in 5 6 7 8; do
    echo -n "$level "
    cat norms-$level
done > norms

if cat <<EOF | gnuplot; then :
set term postscript eps color enhanced lw 2 18 solid

#
# electric field
# 

set output 'efield.eps'
set xlabel 'r'
set ylabel 'E'
set key top right
perm=1.
R0=0.1
rhoinic=0.5
rho(x)= x<R0 ?  0. : 0.5*R0*R0*rhoinic/perm/x
set samples 1000
plot [0:1]\
          'prof-5' t 'R_{o}/h_{min} = 3.2', \
          'prof-6' t 'R_{o}/h_{min} = 6.4', \
          'prof-7' t 'R_{o}/h_{min} = 12.8', \
          'prof-8' t 'R_{o}/h_{min} = 25.6', \
          rho(x) t 'Analytical solution'

#
# error on electric field
#

set output 'error.eps'
set logscale
set xlabel 'R_{o}/h_{min}'
set ylabel 'Error norms on Ey'
set xtics 3.2,2,102.4
set key top right

plot [2.5:30]'norms.ref' u (R0*2**\$1):10 w l t 'Max (ref)', \
             'norms' u (R0*2**\$1):10 w p t 'Max', \
             'norms.ref' u (R0*2**\$1):8 w l t 'L2 (ref)',\
             'norms' u (R0*2**\$1):8 w p t 'L2',\
             .025/x t 'first-order'

#
# error on total amount of charge
# 

set output 'charge.eps'
set logscale
set xlabel 'R_{o}/h_{min}'
set ylabel 'Q (%  Error)'
set xtics 3.2,2,102.4
set key top right

plot [2.5:30]'< tail -q -n 1 rhoe-*' u (R0*2**(5+\$0)):3 w lp t '', 20./x**2 t 'second-order'
EOF
else
    exit 1
fi

# check charge conservation
for level in 6 7 8; do
    if awk 'BEGIN{ charge = 0. }{
      if (charge == 0.)
        charge = $2;
      else if ($2 != charge)
        exit (1);
    }' < rhoe-$level; then :
    else
	exit 1
    fi
done

# check error norms
if cat <<EOF | python ; then :
from check import *
if (Curve('norms',1,8) - Curve('norms.ref',1,8)).max() > 1e-5 or \
   (Curve('norms',1,10) - Curve('norms.ref',1,10)).max() > 1e-5:
    print (Curve('norms',1,8) - Curve('norms.ref',1,8)).max()
    print (Curve('norms',1,10) - Curve('norms.ref',1,10)).max()
    exit(1)
EOF
else
   exit 1
fi
