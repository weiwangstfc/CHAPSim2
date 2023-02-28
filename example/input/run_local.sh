#!/bin/bash
nohup mpirun -np 6 ../../bin/CHAPSim  > OUTPUT_$(date +%Y-%m-%d_%H.%M).log &
