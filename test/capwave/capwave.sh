levels="3 4 5 6"

if test x$donotrun != xtrue; then
    for level in $levels; do
	if sed "s/LEVEL/$level/g" < $1 | gerris2D -; then
	    :
	else
	    exit 1;
	fi
    done
fi

rm -f convergence
for level in $levels; do
    if awk -v level=$level 'BEGIN {s = 0.; n = 0; } {
          t = $1; y = $2;
          getline < "prosperetti"
          s += (y - $2)*(y - $2);
          n += 1;
        }
        END {
          s = sqrt (s/n)/0.01;
          printf ("%d %g\n", level, s);
        }' < wave-$level >> convergence; then
	:
    else
	exit 1;
    fi
done

awk 'BEGIN{first=1}{ 
  if (first) printf("Gerris &\n%.5f",$2);
  else printf(" &\n%.5f",$2);
  first=0;
}' < convergence > convergence.tex

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20
    set output 'amplitude.eps'
    set xlabel 'tau'
    set ylabel 'Relative amplitude'
    plot 'prosperetti' w l t "Prosperetti", 'wave-6' every 10 w p t "Gerris"
EOF
else
    exit 1
fi

if test -f clsvof.tex ; then
    cp convergence.tex gerris.tex
    echo " &\n0.00060" >> gerris.tex
    if cat <<EOF | gnuplot ; then :
      set term postscript eps color lw 3 solid 20
      set output 'convergence.eps'
      set xlabel 'Number of grid points'
      set ylabel 'Relative RMS error'
      set logscale y
      set logscale x 2
      set grid
      plot [5:200][1e-4:1]'gerris.tex' u (2**(\$0 + 2)):1 t "Gerris" w lp, 'prost.tex' u (2**(\$0 + 2)):1 t "PROST" w lp, 'markers.tex' u (2**(\$0 + 2)):1 t "Markers" w lp, 'clsvof.tex' u (2**(\$0 + 2)):1 t "CLSVOF" w lp, 'surfer.tex' u (2**(\$0 + 2)):1 t "Surfer" w lp, 2./x**2 t "Second order"
EOF
    else
	exit 1
    fi
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('convergence',1,2) - Curve('convergence.ref',1,2)).max() > 1e-5:
    exit(1)
EOF
else
   exit 1
fi
