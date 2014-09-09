if test x$donotrun != xtrue; then
    rm -f r0 r1 r2

    for r in 0 1 2; do
	for level in 5 6 7; do
	    if sed "s/LEVEL/$level/g" < $1 | \
		sed "s/BOX/$r/g" | \
		gerris2D - | \
		awk -v level=$level '{
                  print level " " $7 " " $9
                }' >> r$r; then :
	    else
		exit 1
	    fi
	done
    done
fi

if cat <<EOF | python > minion1.tex; then :
from check import *
from sys import *
from math import *

print r"""\begin{tabular}{|c|c|c|c|c|c|}\hline
        & \multicolumn{5}{c|}{\$L_2\$} \\\ \hline
        & \$L=5\$   & \$O_2\$ & \$L=6\$    & \$O_2\$ & \$L=7\$  \\\ \hline"""

def order(r,color='black'):
   for i in range(0,len(r.l)-1):
     y0,y1 = r.l[i][1],r.l[i+1][1]
     print '& {\color{%s}%.2e} & {\color{%s}%4.2f}' % (color, y0, color, log(y0/y1)/log(2.)),
   print '& {\color{%s}%.2e}' % (color, r.l[i+1][1]), r'\\\'

print 'Uniform',
order(Curve('r0',1,2))
order(Curve('r0.ref',1,2),'blue')
print '\$r=1\$',
order(Curve('r1',1,2))
order(Curve('r1.ref',1,2),'blue')
print '\$r=2\$',
order(Curve('r2',1,2))
order(Curve('r2.ref',1,2),'blue')

print r"""\hline
        & \multicolumn{5}{c|}{\$L_\infty\$} \\\ \hline
        & \$L=5\$   & \$O_\infty\$  & \$L=6\$   & \$O_\infty\$  & \$L=7\$ \\\ \hline"""

print 'Uniform',
order(Curve('r0',1,3))
order(Curve('r0.ref',1,3),'blue')
print '\$r=1\$',
order(Curve('r1',1,3))
order(Curve('r1.ref',1,3),'blue')
print '\$r=2\$',
order(Curve('r2',1,3))
order(Curve('r2.ref',1,3),'blue')

print r"\hline\end{tabular}"
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
for r in ['r0','r1','r2']:
  if (Curve(r,1,2) - Curve(r+'.ref',1,2)).max() > 1e-6 or\
     (Curve(r,1,3) - Curve(r+'.ref',1,3)).max() > 1e-6:
    exit(1)
EOF
else
   exit 1
fi
