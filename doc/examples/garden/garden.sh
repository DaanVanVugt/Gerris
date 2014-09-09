#!/bin/sh

if  gerris2D -DMINLEVEL=6 -DNTHETA=24  garden.gfs &&
    gerris2D -DMINLEVEL=0 -DNTHETA=24  garden.gfs &&
    gerris2D -DMINLEVEL=0 -DNTHETA=60  garden.gfs &&
    gerris2D -DMINLEVEL=0 -DNTHETA=120 garden.gfs; then :
else
    exit 1
fi

for i in 6-24 0-24 0-60 0-120; do
    echo "Save end-$i.gnu { format = Gnuplot }" | gfsview-batch2D end-$i.gfs.gz end.gfv
done

for i in 0 24 72 120; do
    echo "Save mesh-$i.gnu { format = Gnuplot }" | gfsview-batch2D sim-0-120-$i.gfs.gz mesh.gfv
done

cat <<EOF | gnuplot
set term postscript eps lw 1 solid 10

set output 'end.eps'
set multiplot
set size 0.5,0.5
set origin 0,0.5
unset key
set xtics 0,1000,4000
set ytics 0,1000,3000
set title 'Non-adaptive 24 directions'
plot [-500:4000][-500:3000]'end-6-24.gnu' u (\$1+2000.):(\$2+2000.) w l
set origin 0.5,0.5
set title 'Adaptive 24 directions'
plot [-500:4000][-500:3000]'end-0-24.gnu' u (\$1+2000.):(\$2+2000.) w l
set origin 0,0
set title 'Adaptive 60 directions'
plot [-500:4000][-500:3000]'end-0-60.gnu' u (\$1+2000.):(\$2+2000.) w l
set origin 0.5,0
set title 'Adaptive 120 directions'
plot [-500:4000][-500:3000]'end-0-120.gnu' u (\$1+2000.):(\$2+2000.) w l
unset multiplot

set output 'mesh.eps'
set size 1,1
set origin 0,0
set multiplot
set size 0.5,0.5
set origin 0,0.5
unset key
set xtics 0,1000,4000
set ytics 0,1000,3000
set title 't = 0'
plot [-500:4000][-500:3000]'mesh-0.gnu' u (\$1+2000.):(\$2+2000.) w l
set origin 0.5,0.5
set title 't = 1 day'
plot [-500:4000][-500:3000]'mesh-24.gnu' u (\$1+2000.):(\$2+2000.) w l
set origin 0,0
set title 't = 3 days'
plot [-500:4000][-500:3000]'mesh-72.gnu' u (\$1+2000.):(\$2+2000.) w l
set origin 0.5,0
set title 't = 5 days'
plot [-500:4000][-500:3000]'mesh-120.gnu' u (\$1+2000.):(\$2+2000.) w l
unset multiplot

EOF

cpu_6_24=`awk '{cpu=$8}END{printf ("%.0f", cpu);}' < log-6-24`
cpu_0_24=`awk '{cpu=$8}END{printf ("%.0f", cpu);}' < log-0-24`
cpu_0_60=`awk '{cpu=$8}END{printf ("%.0f", cpu);}' < log-0-60`
cpu_0_120=`awk '{cpu=$8}END{printf ("%.0f", cpu);}' < log-0-120`

cat <<EOF > cpu.tex
\\begin{tabular}{c|c|c}
Adaptivity & \\# directions & CPU time (seconds)\\\\\\hline
No & 24 & $cpu_6_24 \\\\
Yes & 24 & $cpu_0_24 \\\\
Yes & 60 & $cpu_0_60 \\\\
Yes & 120 & $cpu_0_120
\\end{tabular}
EOF
