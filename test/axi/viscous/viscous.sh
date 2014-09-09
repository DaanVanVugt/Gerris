if test x$donotrun != xtrue; then
    awk 'BEGIN{ for (x = 0.5; x <= 3.; x += 1./256.) print x, 0., 0.;}' > axis
    for Re in 100; do
	if gerris2D -DLEVEL=12 -DRE=$Re viscous.gfs; then :
	else
	    exit 1
	fi
	if gfs2oogl2D -c P -o -i < end-12-$Re.gfs | \
	    awk '{print 3.14159265359 - atan2($2,$1),$4*2.}' | \
	    sort -k 1,2 > cp-12-$Re; then :
	else
	    exit 1
	fi
    done
fi

if echo "Save isolines.eps { format = EPS }" | gfsview-batch2D  end-12-100.gfs isolines.gfv; then :
else
    exit 1
fi

if cat <<EOF | gnuplot ; then :
set term postscript eps color lw 3 solid 20
set key bottom right
set pointsize 1.5
set output 'length.eps'
set xlabel 'Reynolds number'
set ylabel 'Recirculation length'
plot [0:320][0:]'fadlun' t 'Fadlun et al. (2000)', 'fornberg' t 'Fornberg (1988)', 'zhang' t 'Zhang & Zheng (2007)' pt 2, 'blanco-1995' u 1:(\$2/2.) t 'Blanco & Magnaudet (1995)', 'masliyah-1970' u 1:(\$2/2.) t 'Masliyah & Epstein (1970)' pt 8 lt 7, 'Re-12' u 1:3 smooth csplines w l t '' lt 5, 'Re-12' u 1:3 t 'Gerris' lt 5 pt 1

set key top right
set xlabel 'Angle'
set ylabel 'Cp'
set output 'Cp.eps'
plot 'fadlun-cp-100' u (\$1*180./pi):2 w l t 'Fadlun et al., Re = 100', 'cp-12-100' u (\$1*180./pi):2 w l t 'Gerris, Re = 100', 'fadlun-cp-200' u (\$1*180./pi):2 w l t 'Fadlun et al., Re = 200', 'cp-12-200' u (\$1*180./pi):2 w l t 'Gerris, Re = 200'
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('cp-12-100',1,2) - Curve('fadlun-cp-100',1,2)).norm2() > 1.6e-2:
    print (Curve('cp-12-100',1,2) - Curve('fadlun-cp-100',1,2)).norm2()
    exit(1)
EOF
else
   exit 1
fi
