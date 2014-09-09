if test x$donotrun != xtrue; then
    # generate CGD file from SWASHES file
    (
	echo "1 x"
	echo "1000"
	awk '{ if (substr($1,1,1) != "#") printf ("%f ", $1); }' < swashes
	echo ""
	awk '{ if (substr($1,1,1) != "#") print $6; }' < swashes
    ) > href.cgd

    for level in 5 6 7 8; do
	if gerris2D -DLEVEL=$level shock.gfs; then :
	else
	    exit 1
	fi
    done > error
fi

if gnuplot <<EOF; then :
set term postscript eps color lw 2 18 solid enhanced
set output 'error.eps'
set logscale
set xtics 16,2,512
set key spacing 1.5 bottom left
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f1(x)=a1+b1*x
fit f1(x) 'error' u (log(2**\$1)):(log(\$2)) via a1,b1
f2(x)=a2+b2*x
fit f2(x) 'error' u (log(2**\$1)):(log(\$3)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'error' u (log(2**\$1)):(log(\$4)) via am,bm
set xlabel 'Spatial resolution'
set ylabel 'Relative error norms'
plot [22:384][1e-4:]'error' u (2**\$1):4 t 'Max' w lp ps 2, exp(fm(log(x))) t ftitle(am,bm), \
                    'error' u (2**\$1):3 t 'L2' w lp ps 2,  exp(f2(log(x))) t ftitle(a2,b2), \
                    'error' u (2**\$1):2 t 'L1' w lp ps 2,  exp(f1(log(x))) t ftitle(a1,b1) lt 0
EOF
else
    exit 1
fi

if python <<EOF ; then :
from check import *
from sys import *
if (Curve('error',1,3) - Curve('error.ref',1,3)).max() > 1e-5:
    print (Curve('error',1,3) - Curve('error.ref',1,3)).max()
    exit(1)
EOF
else
   exit 1
fi
