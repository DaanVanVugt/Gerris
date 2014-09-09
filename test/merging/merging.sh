if test x$donotrun != xtrue; then
    if sed "s/LEVEL/8/g" < $1 | \
       sed "s/SIM/sim-8/g" | \
       gerris2D - | gfsview-batch2D; then :
    else
	exit 1
    fi

    for level in 6 7 9; do
	if sed "s/LEVEL/$level/g" < $1 | \
           sed "s/SIM/sim-$level/g" | \
           gerris2D - > /dev/null; then :
	else
	    exit 1
	fi
    done

    for level in 6 7 8 9; do
	if sed "s/LEVEL/$level/g" < $1 | \
           sed "s/AdaptVorticity/# AdaptVorticity/g" | \
           sed "s/SIM/simc-$level/g" | \
           gerris2D - > /dev/null; then :
	else
	    exit 1
	fi
    done
fi

for s in sim simc; do
    rm -f $s.err
    for level in 6 7 8; do
	level1=`expr $level + 1`
	echo -n "$level " >> $s.err
	if gfscompare2D -v $s-$level $s-$level1 U 2>&1 | \
	    awk '{if ($1 == "total") print $6 " " $8;}' >> $s.err; then :
	else
	    exit 1
	fi
    done
done

if cat <<EOF | python > convergence.tex ; then :
from check import *
from sys import *
from math import *

print r"""\begin{tabular}{|c|c|c|c|c|c|}\hline
Domain   & \multicolumn{5}{c|}{\$L_2\$}\\\ \hline
         & \$L=6\$   & \$O_2\$ & \$L=7\$    & \$O_2\$ & \$L=8\$  \\\ \hline"""

def order(r,color='black'):
   for i in range(0,len(r.l)-1):
     y0,y1 = r.l[i][1],r.l[i+1][1]
     print '& {\color{%s}%.2e} & {\color{%s}%4.2f}' % (color, y0, color, log(y0/y1)/log(2.)),
   print '& {\color{%s}%.2e}' % (color, r.l[i+1][1]), r'\\\'

print 'Circle',
order(Curve('simc.err',1,2))
order(Curve('simc.err.ref',1,2), 'blue')
print 'Adaptive',
order(Curve('sim.err',1,2))
order(Curve('sim.err.ref',1,2), 'blue')

print r"""\hline
Domain   & \multicolumn{5}{c|}{\$L_\infty\$} \\\ \hline
         &  \$L=6\$   & \$O_\infty\$ & \$L=7\$   & \$O_\infty\$ & \$L=8\$ \\\ \hline"""

print 'Circle',
order(Curve('simc.err',1,3))
order(Curve('simc.err.ref',1,3), 'blue')
print 'Adaptive',
order(Curve('sim.err',1,3))
order(Curve('sim.err.ref',1,3), 'blue')

print r"\hline\end{tabular}"
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('sim.err',1,2) - Curve('sim.err.ref',1,2)).max() > 1e-6 or\
   (Curve('simc.err',1,3) - Curve('simc.err.ref',1,3)).max() > 1e-6:
  exit(1)
EOF
else
   exit 1
fi
