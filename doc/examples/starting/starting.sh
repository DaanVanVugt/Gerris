#!/bin/bash

tail -n +2 $1.dat | shapes - > $1.gts
gerris2D -m -DAEROFOIL=$1 -DINCIDENCE=$2 starting.gfs
