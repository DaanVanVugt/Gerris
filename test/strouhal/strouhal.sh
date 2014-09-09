if test x$donotrun != xtrue; then 
    shapes ellipse > cylinder.gts
    if gerris2D strouhal.gfs > strouhal.res; then :
    else
	exit 1
    fi
fi

if echo "Save stdout { width = 800 height = 200 }" | \
    gfsview-batch2D end.gfs strouhal.gfv | \
    convert -colors 256 ppm:- vort.eps; then :
else
    exit 1
fi

if cat <<EOF | gnuplot ;then :
   set term pos enhanced eps solid color lw 2
   set out 'forces.eps'
   set xl "Time"
   set yl "Cd/Cl"
   set yr [-1.25:2.1]
   set xr [0:112]
   plot 'forces.dat' u (\$1/0.00625):((\$2+\$5)/0.00625) w l t 'Gerris Moving - Cd',\
        'forces.dat' u (\$1/0.00625):((\$3+\$6)/0.00625) w l t 'Gerris Moving - Cl'

EOF
else
   exit 1
fi

if cat <<EOF | gnuplot ;then :
   set term pos enhanced eps solid color lw 2
   set out 'strouhal.eps'
   set xl "Reynolds number"
   set yl "St"
   set yr [0.19:0.24]
   set xr [180:520]
   set key bottom
   plot 'static.ref' u 1:3 w lp t "Gerris Static Low Resolution",\
        'moving.ref' u 1:3 pt 5 t "Gerris Moving Low Resolution",\
        'static.ref' u 1:2 w lp t "Gerris Static High Resolution",\
        'moving.ref' u 1:2 w lp t "Gerris Moving High Resolution",\
        'strouhal.res' pt 5 lc 2 t ""
EOF
else
   exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('strouhal.res',1,2) - Curve('strouhal.ref',1,2)).max() > 0.01 :
    exit(1)
EOF
else
   exit 1
fi
