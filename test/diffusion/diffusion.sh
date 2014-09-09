if gerris2D diffusion.gfs; then
    mv end.gfs end-explicit.gfs
    mv prof prof-explicit
else
    exit 1
fi

if sed -e 's/SourceDiffusionExplicit T/#/' \
       -e 's/#    SourceDiffusion/    SourceDiffusion/' < diffusion.gfs | \
   gerris2D -; then
    mv end.gfs end-implicit.gfs
    mv prof prof-implicit
else
    exit 1
fi

if cat <<EOF | gnuplot; then :
    set term postscript eps color lw 3 solid 20
    set output 'profile.eps'
    set xlabel 'r'
    set ylabel 'T'
    set key top left
    plot 'prof-explicit' t 'explicit', 'prof-implicit' t 'implicit'
EOF
else
    exit 1
fi

if gfscompare2D -v end-explicit.gfs end-implicit.gfs T 2>&1 | awk '{
  if ($1 == "total" && ($6 > 3e-3 || $8 > 4e-3)) {
    print $0
    exit (1)
  }
}'; then :
else
    exit 1
fi
