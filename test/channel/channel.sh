if test x$donotrun != xtrue; then
    shapes channel | transform --revert --scale 4 --tx 1.5 > channel.gts
    for level in 5 6 7; do
	if sed "s/LEVEL/$level/g" < $1 | \
           gerris2D -; then :
	else
	    exit 1
	fi
    done
fi

for v in U V; do
    rm -f order$v orderf$v
    for level in 5 6; do
	level1=`expr $level + 1`
	echo -n "$level " >> order$v
	if gfscompare2D -v sim-$level sim-$level1 $v 2>&1 | \
	    awk '{if ($1 == "total") print $4 " " $6 " " $8;}' >> order$v; then :
	else
	    exit 1
	fi
	echo -n "$level " >> orderf$v
	if gfscompare2D -f 5 -v sim-$level sim-$level1 $v 2>&1 | \
	    awk '{if ($1 == "total") print $4 " " $6 " " $8;}' >> orderf$v; then :
	else
	    exit 1
	fi
    done
done

if cat <<EOF | python > convergence.tex; then :
from check import *
from sys import *
from math import *

for component,variable in [('x','U'),('y','V')]:
  print r"""\begin{table}[htbp]
  \caption{"""
  print r"\label{channel-" + component + "}"
  print r"Errors and convergence rates for the \$"+component+r"\$-component of the velocity.}"
  print r"""\begin{center}
  \begin{tabular}{||l|c|c|c||c|c|c||} \hline
           & \multicolumn{3}{c||}{All cells} & \multicolumn{3}{c||}{Full 128 cells} \\\ \hline
           & 128-256  & Rate & 256-512  & 128-256  & Rate & 256-512  \\\ \hline"""

  for i,name in [(2,r"\$L_1\$"),(3,r"\$L_2\$"),(4,r"\$L_\infty\$")]:
    a=Curve('order'+variable,1,i)
    b=Curve('orderf'+variable,1,i)
    print name,
    print "& %.2e & %4.2f & %.2e & %.2e & %4.2f & %.2e" % (\
    a.l[0][1], log(a.l[0][1]/a.l[1][1])/log(2.), a.l[1][1], \
    b.l[0][1], log(b.l[0][1]/b.l[1][1])/log(2.), b.l[1][1]),
    print r"\\\"

    a=Curve('order'+variable+'.ref',1,i)
    b=Curve('orderf'+variable+'.ref',1,i)
    print "& {\color{blue}%.2e} & {\color{blue}%4.2f} & {\color{blue}%.2e} & {\color{blue}%.2e} & {\color{blue}%4.2f} & {\color{blue}%.2e}" % (\
    a.l[0][1], log(a.l[0][1]/a.l[1][1])/log(2.), a.l[1][1], \
    b.l[0][1], log(b.l[0][1]/b.l[1][1])/log(2.), b.l[1][1]),
    print r"\\\"

  print r"\hline"
  print r"""\end{tabular}
  \end{center}
  \end{table}"""
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *

for f in ['orderU','orderV','orderfU','orderfV']:
   if (Curve(f,1,2) - Curve(f+'.ref',1,2)).max() > 1e-6 or\
      (Curve(f,1,3) - Curve(f+'.ref',1,3)).max() > 1e-6 or\
      (Curve(f,1,4) - Curve(f+'.ref',1,4)).max() > 1e-6:
      exit(1)
EOF
else
   exit 1
fi
