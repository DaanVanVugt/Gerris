if gerris2D -m concentration.gfs; then
    mv end.gfs end-explicit.gfs
    mv prof prof-explicit
else
    exit 1
fi

if sed -e 's/SourceDiffusionExplicit C/#/' \
       -e 's/#    SourceDiffusion/    SourceDiffusion/' < concentration.gfs | \
   gerris2D -m -; then
    mv end.gfs end-implicit.gfs
    mv prof prof-implicit
else
    exit 1
fi

if cat <<EOF | gnuplot; then :
    set term postscript eps color lw 3 solid 20
    set output 'profile.eps'
    set xlabel 'r'
    set ylabel 'T1,C'
    set key top left
    plot 'prof-explicit' u 1:2 t 'Tracer', \
         'prof-explicit' u 1:3 t 'C explicit', \
         'prof-implicit' u 1:3 t 'C implicit'
EOF
else
    exit 1
fi

if gfscompare2D -v end-explicit.gfs end-implicit.gfs C 2>&1 | awk '{
  if ($1 == "total" && ($6 > 7.07e-4 || $8 > 2.27e-3)) {
    print $0
    exit (1)
  }
}'; then :
else
    exit 1
fi
