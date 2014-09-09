#!/bin/bash

set -e

disable_hypre() {
	## comment out hypre
	perl -wpl -i -e 's/^(\s*)(GModule hypre)/$1#$2/' $1
}

enable_hypre() {
	## uncomment hypre (if commented)
	perl -wpl -i -e 's/#(GModule hypre)/$1/' $1
}

NP=4
CASEFILE=${1}

[ -f "${CASEFILE}" ] || {
	echo "${CASEFILE} does not exist"
	exit 1
}

## if the hypre module is unavailable
## is there a way to automatically check whether hypre is available?
#disable_hypre ${CASEFILE}
## otherwise
enable_hypre ${CASEFILE}

SPLITFILE="${CASEFILE%.gfs}-split.gfs"
gerris2D --split=3 -m ${CASEFILE} | gerris2D --partition=2 - > ${SPLITFILE}
mpirun -np $NP gerris2D ${SPLITFILE} | gfsview-batch2D bubble.gfv
rm ${SPLITFILE}
