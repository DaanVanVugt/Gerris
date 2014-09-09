#!/bin/bash
export GERRIS=gerris2D
export FILENAME=viscmix

##  RUN JOB
# echo && echo start job on `hostname` at `date` && echo && echo

# split the domain twice to get enough boxes to redistribute
mpirun ${GERRIS} -s 2 ${FILENAME}.gfs > ${FILENAME}_split.gfs && \

# create the initial partition into 2^2=4 subdomains
mpirun ${GERRIS} -p 2 ${FILENAME}_split.gfs > ${FILENAME}_partitioned.gfs && \

# run the parallel simulation on 4 processors
mpirun -np 4 ${GERRIS} ${FILENAME}_partitioned.gfs

# echo && echo && echo job finished on `hostname` at `date`.
