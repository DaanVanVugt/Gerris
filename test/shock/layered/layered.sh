if test x$donotrun != xtrue; then
    for l in 1 2 5 15; do
	if gerris2D -DLEVEL=8 -DNL=$l layered.gfs; then :
	else
	    exit 1
	fi
    done
fi

if gnuplot <<EOF; then :
set term postscript eps lw 2 18 color enhanced solid

set output 'prof.eps'
set xlabel 'x'
set ylabel 'z'
set xrange [0:21]
set yrange [0:1]
set xtics 0,2,20
plot 0.2*(1. - 1./(5.75/2.)**2*(x - 10.)**2) t 'topography' w l, \
     'prof-8-1' w l t '1 layer', \
     'prof-8-2' w l t '2 layers', \
     'prof-8-5' w l t '5 layers', \
     'prof-8-15' w l t '15 layers'

set term postscript eps lw 2 24 color enhanced solid

set xrange [0:0.6]
set yrange [0.5:2.5]
set key bottom right
set pointsize 2
set xtics auto

set output 'uprof-10.eps'
set xlabel 'z - z_b'
set ylabel 'u'
plot 'uprof-10-8-2' t '2 layers', 'uprof-10-8-5' t '5 layers', 'uprof-10-8-15' t '15 layers'

set output 'uprof-15.eps'
plot 'uprof-15-8-2' t '2 layers', 'uprof-15-8-5' t '5 layers', 'uprof-15-8-15' t '15 layers'

set output 'uprof-20.eps'
plot 'uprof-20-8-2' t '2 layers', 'uprof-20-8-5' t '5 layers', 'uprof-20-8-15' t '15 layers'

set term postscript eps color enhanced solid lw 1 16
set output 'u.eps'
set pm3d
set pm3d map interpolate 10,1
set xlabel 'x'
set ylabel 'z'
unset key
# jet colormap
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,\
     0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0,\
     0.875 0.9333 0 0, 1 0.498 0 0 )
splot [0:21][0:1]'< awk -v name="U" -f field.awk end-8-15' u 1:2:3
EOF
else
    exit 1
fi
fixbb u.eps

if python <<EOF; then :
from check import *
from sys import *
if (Curve('prof-8-15',1,2) - Curve('prof.ref',1,2)).max() > 1e-3:
    print (Curve('prof-8-15',1,2) - Curve('prof.ref',1,2)).max()
    exit(1)
EOF
else
   exit 1
fi
