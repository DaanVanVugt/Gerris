if test x$donotrun != xtrue; then
    gerris2D axiadvection.gfs | gfsview-batch2D -s vof.gfv > vof.gnu
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20
    set output 'vof.eps'
    set xlabel 'z'
    set ylabel 'r'
    unset key
    set xtics -0.5,.25,0.5
    set size ratio -1
    plot [-0.5:0.5][0:1]'vof.gnu' u 1:2 w l, 'vectors.gnu' u 1:2 w l
EOF
else
    exit 1
fi

if awk '
BEGIN { min = 1000.; max = -1000.; }{ 
  if ($5 < min) min = $5; 
  if ($5 > max) max = $5; 
}
END {
  e = 2.*(max - min)/(max + min);
  print "VOF:", e;
  if (e > 4e-4)
    exit (1);
}' < srt; then :
else
    exit 1
fi

if awk '
BEGIN { min = 1000.; max = -1000.; }{ 
  if ($5 < min) min = $5; 
  if ($5 > max) max = $5; 
}
END {
  e = 2.*(max - min)/(max + min);
  print "Standard:", e;
  if (e > 4e-7)
    exit (1);
}' < srt1; then :
else
    exit 1
fi
