#!/bin/bash
nohup mpirun -np 4 ../../bin/CHAPSim <input_chapsim.ini > OUTPUT_$(date +%Y-%m-%d_%H.%M).log &
