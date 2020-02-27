#!/bin/bash

cflags='-O3'

gfortran -c overlapFINAL.F90
gfortran -c kinetic_monte_carlo_snapshots.f90 
gfortran -c snapshot_KMC.f90
gfortran $cflags overlapFINAL.o kinetic_monte_carlo_snapshots.o snapshot_KMC.o -o snapshot_KMC.x -L/usr/lib64 -lblas

rm *.o *.mod

