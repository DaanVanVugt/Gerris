# split three times and partition on 4 processors
gerris3D -s 3 atomisation.gfs | gerris3D -p2 - > atomisation-s3-p2.gfs

# run the parallel simulation
mpirun -np 4 gerris3D atomisation-s3-p2.gfs < /dev/null
