if test x$donotrun != xtrue; then
    for model in 0 1 2 3; do
	if sed "s/MODEL/$model/g" < $1 | gerris2D -; then :
	else
	    exit 1
	fi
    done
fi

if cat <<EOF | gnuplot ; then :
set term postscript eps color lw 3 solid 20
set output 'prof.eps'
set xlabel 'r'
set ylabel 'Tangential velocity'
powerlaw(r,N)=r*((0.5/r)**(2./N) - 1.)/((0.5/0.25)**(2./N) - 1.)
hb(r,Rl)=(r > Rl ? 0. : r*sqrt(2.)*0.12*0.12/(4.*0.0672*0.0672)*(3./4.+(Rl/r)**4/4.-(Rl/r)**2+log(Rl/r)))
bingham(r,Rl)=(r > Rl ? 0. : r*sqrt(2.)*10./4.*((Rl/r)**2-2.*log(Rl/r)-1.))
plot [0.25:0.5][0:0.25]powerlaw(x,1.) t "Newtonian", 'prof-0' w p ps 2 pt 9 t "",\
               powerlaw(x,0.5) t "Power law", 'prof-1' w p ps 2 pt 9 t "",\
               hb(x,0.4637) t "Herschel-Bulkley", 'prof-2' w p ps 2 pt 9 t "",\
               bingham(x,0.34924) t "Bingham", 'prof-3' w p ps 2 pt 9 t ""
EOF
else
   exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
print (Curve('prof-0',1,2) - Curve('prof-0.ref',1,2)).norm2(),\
   (Curve('prof-1',1,2) - Curve('prof-1.ref',1,2)).norm2(),\
   (Curve('prof-2',1,2) - Curve('prof-2.ref',1,2)).norm2(),\
   (Curve('prof-3',1,2) - Curve('prof-3.ref',1,2)).norm2()
if (Curve('prof-0',1,2) - Curve('prof-0.ref',1,2)).norm2() > 3.6e-4 or \
   (Curve('prof-1',1,2) - Curve('prof-1.ref',1,2)).norm2() > 6.3e-4 or \
   (Curve('prof-2',1,2) - Curve('prof-2.ref',1,2)).norm2() > 21e-4 or \
   (Curve('prof-3',1,2) - Curve('prof-3.ref',1,2)).norm2() > 22e-4:
    exit(1)
EOF
else
   exit 1
fi
