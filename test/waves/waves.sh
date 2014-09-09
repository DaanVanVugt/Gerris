if test x$donotrun != xtrue; then
    rm -f correlation res-* sim-*
    for level in 5 6 7; do
	if sed "s/LEVEL/$level/g" < $1 | gerris2D - | \
	    awk '{ print $1 " " $5 }' > res-$level && \
	    awk -v l=$level 'BEGIN { min1 = 0. } {
         if ($2 > min1) {
           theta = $1;
           min1 = $2;
         }
       } END {
           printf ("%d\t\t%g\t\t\t%g\n", l, theta, min1);
       }' < res-$level >> correlation; then :
	else
	    exit 1
	fi
    done
fi

echo "Save solution.eps { format = EPS }" | gfsview-batch2D sim-7 solution.gfv

awk 'BEGIN {
  print "\\begin{tabular}{|c|c|c|}"
  print "\\hline Level & Maximum $C$ & Angle of max $C$ \\\\ \\hline"
}{
  print $1 " & " $3 " & " $2 " \\\\"
}END {
  print "\\hline"
  print "\\end{tabular}"
}' < correlation > correlation.tex

if cat <<EOF | python ; then :
from check import *
from sys import *
if Curve('correlation',1,3).max() > 10.:
    exit(1)
if (Curve('correlation',1,2) - Curve('correlation.ref',1,2)).max() > 0. or\
   (Curve('correlation.ref',1,3) - Curve('correlation',1,3)).max() > 0.:
    exit(1)
EOF
else
   exit 1
fi
