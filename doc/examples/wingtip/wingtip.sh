#!/bin/sh

tail -n +2 $1.dat | shapes - > $1.gts
gerris3D -m -DAEROFOIL=$1 -DINCIDENCE=$2 wingtip.gfs
