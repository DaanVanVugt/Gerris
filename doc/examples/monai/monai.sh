# get bathymetry data
if test ! -f Benchmark_2_Bathymetry.txt; then
    wget http://isec.nacse.org/workshop/2004_cornell/data/Benchmark_2_Bathymetry.txt
fi

# number of lines in the file (ignoring the first line which is a header)
np=`awk 'FNR>1 && NF == 3' Benchmark_2_Bathymetry.txt | wc -l`

# triangulate the data points
(echo "$np 0 0" && awk 'FNR>1 && NF == 3 { print $1,$2,-$3 }' Benchmark_2_Bathymetry.txt) | delaunay -r -v > bathy.gts

# get input data
if test ! -f Benchmark_2_input.txt; then
    wget http://isec.nacse.org/workshop/2004_cornell/data/Benchmark_2_input.txt
fi

# number of lines in the file (ignoring the first line which is a header)
np=`awk 'FNR>1 && NF == 2' Benchmark_2_input.txt | wc -l`

# create CGD input file
cat <<EOF > input.cgd
1 t
$np
EOF
awk 'FNR>1 && NF == 2 { printf ("%s ",$1) }' Benchmark_2_input.txt >> input.cgd
awk 'FNR>1 && NF == 2 { print $2 }' Benchmark_2_input.txt >> input.cgd

# run the simulation
if xdpyinfo > /dev/null 2>&1; then 
    GFSVIEW=gfsview2D
else 
    GFSVIEW=gfsview-batch2D
fi

if gerris2D -m monai.gfs | $GFSVIEW 3D.gfv; then :
else
    exit 1
fi

# generate graphics
for i in 10 12 14 16 18 20; do
    if test $i = 18; then
	probe=probe.gfv
    else
	probe=""
    fi
    echo "Save stdout { width = 602 height = 582 }" | \
	gfsview-batch2D sim-$i.gfs.gz leveque.gfv $probe | convert ppm:- eps2:fig4-$i.eps
    echo "Save stdout { width = 602 height = 582 }" | \
	gfsview-batch2D sim-$i.gfs.gz mesh.gfv | convert ppm:- eps2:mesh-$i.eps
done
echo "Save stdout { width = 1280 height = 960 }" | \
    gfsview-batch2D sim-18.gfs.gz 3D.gfv | convert ppm:- eps2:monai.eps

# get experimental probe data
if test ! -f output_ch5-7-9.xls; then
    wget http://isec.nacse.org/workshop/2004_cornell/data/benchmark2/output_ch5-7-9.xls
fi

# convert excel crap to plain text (requires catdoc, install with 'sudo apt-get install catdoc')
xls2csv -c' ' output_ch5-7-9.xls | sed 's/"//g' | awk 'FNR>1' > output_ch5-7-9.txt

gnuplot <<EOF
set term postscript eps color lw 2 20 solid
set xlabel "Time (s)"
set ylabel "Elevation (cm)" 
set key top left
set output 'p5.eps'
plot [0:22.5]'output_ch5-7-9.txt' u 1:2 pt 6 ps 0.5 t 'Experiment', 'p5' u 1:((\$9)*100.) w l t 'Gerris'
set output 'p7.eps'
plot [0:22.5]'output_ch5-7-9.txt' u 1:3 pt 6 ps 0.5 t 'Experiment', 'p7' u 1:((\$9)*100.) w l t 'Gerris'
set output 'p9.eps'
plot [0:22.5]'output_ch5-7-9.txt' u 1:4 pt 6 ps 0.5 t 'Experiment', 'p9' u 1:((\$9)*100.) w l t 'Gerris'
EOF

# generate overhead comparison movie
# this link is broken
# wget http://isec.nacse.org/workshop/2004_cornell/data/benchmark2/overhead.avi
# use this one instead
if test ! -f experiment.mpg; then
    wget -O experiment.mpg http://www.amath.washington.edu/~rjl/catalina04/overhead.mpg
fi

# extract individual frames and change contrast of experimental movie
ffmpeg -i experiment.mpg frame-%03d.png
mogrify -modulate 300,100,100 -contrast -rotate 90 -chop 0x10 -geometry x600 frame-*.png

# extract individual frames of simulation
ffmpeg -i overhead.mpg sframe-%03d.png

# combine frames
for f in frame-*.png; do
    montage -geometry +0+0 $f s$f -resize 400x600! -depth 8 png:m$f
    echo -n -e '\r'$f
done
echo ""
ffmpeg -r 10 -f image2 -i mframe-%03d.png -b 1800K comparison.mp4

convert mframe-050.png eps2:comparison.eps

# cleanup
rm -f frame-*.png sframe-*.png mframe-*.png
