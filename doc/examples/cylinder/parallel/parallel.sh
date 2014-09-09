# !/bin/sh

# split the domain twice to get enough boxes to redistribute
gerris2D -s 2 parallel.gfs > parallel-s2.gfs

# create the initial partition into 2^2=4 subdomains
gerris2D -p 2 parallel-s2.gfs > parallel-p2.gfs

# run the parallel simulation on 4 processors, pipe the output to
# gfsview and ppm2mpeg to generate the pid movie
mpirun -np 4 gerris2D parallel-p2.gfs
