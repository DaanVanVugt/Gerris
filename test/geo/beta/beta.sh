if test x$donotrun != xtrue; then
    if gerris3D $1; then :
    else
	exit 1
    fi
    mv -f energy energy-nonlinear
    if sed 's/Refine 6/Refine 6\nAdvectionParams {scheme = none}/' < $1 |\
       gerris3D -; then :
    else
	exit 1
    fi
fi

if cat <<EOF | gnuplot; then :
    set term postscript eps lw 3 color solid 20
    set output 'energy.eps'
    set xlabel 'Time (days)'
    set ylabel 'Normalised total energy'
    set key bottom left
    plot 'c' t 'C-grid' w lp, 'energy' t 'Gerris (linear)' w lp, 'energy-nonlinear' t 'Gerris' w lp, 'dlw' t 'Delumped LW' w lp, 'lls' t 'LLS' w lp, 'pzm' t 'PZM' w lp, 'llw' t 'Lumped LW' w lp
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if Curve('energy',1,2).max() > 1. or \
   Curve('energy-nonlinear',1,2).max() > 1.:
    exit(1)
if (Curve('energy.ref',1,2) - Curve('energy',1,2)).max() > 1e-2 or \
   (Curve('energy-nonlinear.ref',1,2) - Curve('energy-nonlinear',1,2)).max() > 1e-2:
    exit(1)
EOF
else
   exit 1
fi
