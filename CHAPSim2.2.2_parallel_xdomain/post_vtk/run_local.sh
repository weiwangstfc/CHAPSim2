#!/bin/bash
nohup mpirun -np 1 ./src/post_vtk  > post_$(date +%Y-%m-%d_%H.%M).log &
