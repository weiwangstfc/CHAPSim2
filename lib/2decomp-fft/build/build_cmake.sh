#!/bin/sh
FC=mpif90 
cmake -S ../ -B ./
cmake --build ./
cmake --install ./
