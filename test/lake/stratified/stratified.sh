ratio=10
if test x$donotrun != xtrue; then
    for nl in 8 16 32 64; do
	if gerris2D -DRATIO=$ratio -DNL=$nl -DRE=100 -DLEVEL=5 stratified.gfs; then :
	else
	    exit 1
	fi
    done
fi

for nl in 8 16 32 64; do
    awk -v nl=$nl '{
      z = $3
      u = (z > 3./4. ? (30.*z*z - 47.*z + 18.)/13. : (-2.*z*z + z)/13.)
      e = $4 - u
      if (e < 0.) e = -e;
      if (e > emax) emax = e;
    } END {
      print nl,emax
    }' < uprof-$nl
    awk -f field.awk -v name="DRHO" end-$nl.txt | awk -f thermo.awk > thermo-$nl.txt
done > error
awk -f field.awk -v name="DRHO" end-64.txt > rho-64.txt

if gnuplot <<EOF; then :
set term postscript eps color lw 2 18 solid enhanced
set output 'uprof.eps'
set key top left
set xlabel 'z'
set ylabel 'u'
uprof(z)=(z > 3./4. ? (30.*z*z - 47.*z + 18.)/13. : (-2.*z*z + z)/13.)
plot [0:1] uprof(x) t 'analytical', \
          'uprof-8' u 3:4 t '8 layers', \
          'uprof-16' u 3:4 t '16 layers', \
          'uprof-32' u 3:4 t '32 layers', \
          'uprof-64' u 3:4 t '64 layers'

set output 'free.eps'
set key top left
set xlabel 'x'
set ylabel 'z'
plot 60./13.*1e-4*x+1 t 'analytical', \
     'end-8.txt' u 1:4 t '8 layers', \
     'end-16.txt' u 1:4 t '16 layers', \
     'end-32.txt' u 1:4 t '32 layers', \
     'end-64.txt' u 1:4 t '64 layers'

set output 'thermo.eps'
set key top right
plot -64./13.*1e-3*x+0.75 t 'analytical', \
     'thermo-8.txt' t '8 layers', \
     'thermo-16.txt' t '16 layers', \
     'thermo-32.txt' t '32 layers', \
     'thermo-64.txt' t '64 layers'

set output 'error.eps'
fm(x)=am+bm*x
fit fm(x) 'error' u (log(\$1)):(log(\$2)) via am,bm
set key spacing 1.5 top right
ftitle(a,b) = sprintf("%.2g/x^{%4.2f}", exp(a), -b)
set xlabel 'Spatial resolution'
set ylabel 'Maximum error'
set logscale
set xtics 4,2,128
plot [4:128]'error' u 1:2 t '' w lp ps 2, exp(fm(log(x))) t ftitle(am,bm)

set term postscript eps color enhanced solid lw 1 16
set output 'drho.eps'
reset
set pm3d
set pm3d map interpolate 4,4
set xlabel 'x'
set ylabel 'z'
unset key
# jet colormap
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,\
     0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0,\
     0.875 0.9333 0 0, 1 0.498 0 0 )
splot [-5:5][0:1.01]'rho-64.txt' u 1:2:3
EOF
else
    exit 1
fi

if python <<EOF ; then :
from check import *
from sys import *
if (Curve('error',1,2) - Curve('error.ref',1,2)).max() > 1e-6:
    print (Curve('error',1,2) - Curve('error.ref',1,2)).max()
    exit(1)
EOF
else
   exit 1
fi
