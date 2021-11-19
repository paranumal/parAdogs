#!/bin/bash

rm -f *.vtu

mpirun -np $1 ./paradogsMain $2

paraview --script=paraview_vis.py
