#!/bin/bash
nohup mpirun -np 1 ../../bin/CHAPSim  > OUTPUT_$(date +%Y-%m-%d_%H.%M).log &
