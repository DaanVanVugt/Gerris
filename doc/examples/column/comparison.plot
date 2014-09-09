r0=0.677/6.26
unset xtics
unset ytics
unset border
set size ratio -1
unset key

set term postscript eps color lw 2 enhanced size 5,3.5
set xrange [-5:5]
set yrange [0:7]

set output 'comparison-0.eps'
plot "< awk '{if ($4 > 0) print $0;}' < grains-0.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "#8f8f8f" fs solid 1 noborder, "< awk '{if ($4 == 0) print $0;}' < grains-0.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "black" fs solid 1 noborder, 'snapshot-0.gnu' w l lw 2 lt 1

set output 'comparison-0.66.eps'
plot "< awk '{if ($4 > 0) print $0;}' < grains-0.66.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "#8f8f8f" fs solid 1 noborder, "< awk '{if ($4 == 0) print $0;}' < grains-0.66.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "black" fs solid 1 noborder, 'snapshot-1.65132.gnu' w l lw 2 lt 1

set term postscript eps color lw 2 enhanced size 5,1.75
set xrange [-5:5]
set yrange [0:3.5]

set output 'comparison-0.95.eps'
plot "< awk '{if ($4 > 0) print $0;}' < grains-0.95.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "#8f8f8f" fs solid 1 noborder, "< awk '{if ($4 == 0) print $0;}' < grains-0.95.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "black" fs solid 1 noborder, 'snapshot-2.3769.gnu' w l lw 2 lt 1

set output 'comparison-1.24.eps'
plot "< awk '{if ($4 > 0) print $0;}' < grains-1.24.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "#8f8f8f" fs solid 1 noborder, "< awk '{if ($4 == 0) print $0;}' < grains-1.24.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "black" fs solid 1 noborder, 'snapshot-3.10248.gnu' w l lw 2 lt 1

set term postscript eps color lw 2 enhanced size 10,1
set xrange [-10:10]
set yrange [0:2]

set output 'comparison-1.52.eps'
plot "< awk '{if ($4 > 0) print $0;}' < grains-1.52.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "#8f8f8f" fs solid 1 noborder, "< awk '{if ($4 == 0) print $0;}' < grains-1.52.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "black" fs solid 1 noborder, 'snapshot-3.80304.gnu' w l lw 2 lt 1

set output 'comparison-2.28.eps'
plot "< awk '{if ($4 > 0) print $0;}' < grains-2.28.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "#8f8f8f" fs solid 1 noborder, "< awk '{if ($4 == 0) print $0;}' < grains-2.28.dat" u ($1/r0):($2/r0):($3/r0) w circles lc rgb "black" fs solid 1 noborder, 'snapshot-5.70456.gnu' w l lw 2 lt 1
