if test x$donotrun != xtrue; then
    for order in 1 2 ; do
	if gerris2D -DORDER=$order hexagon.gfs; then :
	else
	    exit 1
	fi
    done

fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20
    set output 'error.eps'
    set xlabel 'Time'
    set ylabel 'Norm2 of error on the velocity field'
    set yr [1e-10:0.1]
    set log y
    plot 'momentumerror-1' t 'first order' w lp, \
         'momentumerror-2' t 'second order' w lp
EOF
else
    exit 1
fi

if echo "Save end-2.eps { format = EPS }" | gfsview-batch2D end-2.gfs hexagon.gfv; then :
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if Curve('tracersum-1',1,2).normi() > 1e-6 or \
   Curve('tracersum-2',1,2).normi() > 5e-6 or \
   Curve('momentumerror-1',1,3).max() > 11e-3 or \
   Curve('momentumerror-2',1,3).max() > 2e-3:
    exit(1)
EOF
else
    exit 1
fi
