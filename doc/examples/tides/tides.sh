#
# Bathymetry
#

# happrox is distributed with GTS. Do 'happrox -h' for more info.
# 'bathymetry' contains the (negative) depth as a function of space i.e. three columns:
# lon (deg) lat (deg) depth (m).
# 0.05 is the relative depth error allowed.
happrox -f -r 1 -c 0.05 < bathymetry | transform --revert > bath.gts

#
# M2 tidal coefficients
#

# 'coefficients' contains the amplitude and phase of the M2
# tide as a function of space i.e. four columns:
# lon (deg) lat (deg) amplitude (m) phase (degree)

# just counts the number of lines in the input file
lines=`wc -l coefficients | awk '{print $1}'`

# Gerris defines the tidal M2 mode as:
# AM2*cos (omega*t) + BM2*sin (omega*t)
#
# We first compute AM2 from the amplitude and phase.
awk -v lines=$lines '
BEGIN {
  print lines " 0 0"
} {
  print $1 " " $2 " " $3*cos($4*3.14159265357/180.)
}' < coefficients | delaunay > AM2.gts

# Now compute BM2
awk -v lines=$lines '
BEGIN {
  print lines " 0 0"
} {
  print $1 " " $2 " " $3*sin($4*3.14159265357/180.)
}' < coefficients | delaunay > BM2.gts

# Run the simulation
gerris3D -m tides.gfs | gfsview3D tides.gfv

# Use batch mode of gfsview to generate figures
echo "Save amplitude.eps { format = EPS }" | gfsview-batch3D end.gfs amplitude.gfv
echo "Save phase.eps { format = EPS }" | gfsview-batch3D end.gfs phase.gfv
echo "Save ellipses.eps { format = EPS }" | gfsview-batch3D end.gfs ellipses.gfv
echo "Save residual.eps { format = EPS }" | gfsview-batch3D end.gfs residual.gfv
