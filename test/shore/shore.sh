levelmin=11
levelmax=14

if test x$donotrun != xtrue; then
    for i in `seq $levelmin 1 $levelmax`; do
	if gerris2D -DLEVEL=$i shore.gfs 2> log-$i; then :
	else
	    exit 1
	fi
    done
fi

if python <<EOF > convergence; then :
from check import *
from sys import *

for i in range($levelmin,$levelmax + 1):
  ref = Curve('t220.csv',1,2)
  sol = Curve('sim-' + str(i) + '-220.txt',1,8)
  diff = ref - sol
  print i,diff.norm1(),diff.normi()

EOF
else
   exit 1
fi

if gnuplot <<EOF ; then :
set term postscript eps color enhanced lw 2 18 solid

set output 'profile.eps'
set xlabel 'x (m)'
set ylabel 'y (m)'
plot [-200:1000][-25:20] \
     -x/10 lc 0 t 'beach', \
     't160.csv' w l t 't = 160 sec', \
     'sim-$levelmax-160.txt' u 1:8 every 6 t '', \
     't175.csv' w l t 't = 175 sec', \
     'sim-$levelmax-175.txt' u 1:8 every 6 t '', \
     't220.csv' w l lc 1 t 't = 220 sec', \
     'sim-$levelmax-220.txt' u 1:8 every 6 t ''

set output 'order.eps'
set logscale
set key spacing 1.5 bottom right
ftitle(a,b) = sprintf("%.2g x^{%4.2f}", exp(a), b)
f2(x)=a2+b2*x
fit f2(x) 'convergence' u (log(60000./2**\$1)):(log(\$2)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'convergence' u (log(60000./2**\$1)):(log(\$3)) via am,bm
set xlabel 'Spatial resolution (m)'
set ylabel 'Error norms (m)'
set xtics 5,2,40
set grid
plot [2:40][0.05:5]\
     'convergence' u (60000./2**\$1):2 t 'L1' w lp ps 2,  exp(f2(log(x))) t ftitle(a2,b2), \
     'convergence' u (60000./2**\$1):3 t 'Max' w lp ps 2, exp(fm(log(x))) t ftitle(am,bm)
EOF
else
    exit 1
fi

if python <<EOF ; then :
from check import *
from sys import *
if (Curve('convergence',1,2) - Curve('convergence.ref',1,2)).max() > 0. or\
   (Curve('convergence',1,3) - Curve('convergence.ref',1,3)).max() > 0.:
    exit(1)
EOF
else
   exit 1
fi