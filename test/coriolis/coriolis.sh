if test x$donotrun != xtrue; then
    if gerris3D -m $1; then :
    else
	exit 1
    fi
fi

if gnuplot <<EOF ; then :
set term pos enhanced eps color solid lw 2 18
set out 'velocity.eps'

set bmargin 4
set rmargin 3.5

set xl 'Time (s)'
set yl 'Velocity Components (m/s)'

set key spacing 3

U=2
V=3
W=1
Phi=3.14159/4.
Omega=7.292e-5
A=U
B=(sin(Phi)*V-cos(Phi)*W)
C=(cos(Phi)*(cos(Phi)*V+sin(Phi)*W))
D=(sin(Phi)*(cos(Phi)*V+sin(Phi)*W))
solu(t)=(A*cos(2*Omega*t)+B*sin(2*Omega*t))
solv(t)=(-A*sin(Phi)*sin(2*Omega*t) + B*sin(Phi)*cos(2*Omega*t) + C)
solw(t)=(A*cos(Phi)*sin(2*Omega*t) - B*cos(Phi)*cos(2*Omega*t) + D)

plot solu(x) t '', solv(x) t '', solw(x) t '',\
     'u.dat' w p lc 1 ps 1.8 t 'U',\
     'v.dat' w p lc 2 ps 1.8 t 'V',\
     'w.dat' w p lc 3 ps 1.8 t 'W'

EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('error.dat',1,2)).max() > 1e-6:
    exit(1)
EOF
else
   exit 1
fi
