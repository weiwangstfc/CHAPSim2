#!/bin/bash

# Define the relative paths to the directories containing the Makefiles
REL_PATH_A="./lib/2decomp_fft_updated/src"
REL_PATH_B="./build"

# Get the absolute paths from the relative paths
PATH_A=$(realpath "$REL_PATH_A")
PATH_B=$(realpath "$REL_PATH_B")

# Ask the user if they want to run 'make clean' in path_B
echo "Do you want to run 'make clean' in CHAPSim $PATH_B? (yes/no): "
read CLEAN_B

# Ask the user for the make target
echo "Enter the make target (e.g., 'make' or 'make all'): "
read MAKE_TARGET

# Check if the user provided a valid target (e.g., 'make' or 'make all')
if [[ -z "$MAKE_TARGET" ]]; then
  echo "No make target provided. Exiting."
  exit 1
fi

# Change to the first directory (A) and run the Makefile
echo "Running Makefile for libs..."
cd "$PATH_A" || { echo "Failed to change directory to $PATH_A"; exit 1; }
make || { echo "Make failed in $PATH_A"; exit 1; }

# If the user wants to run 'make clean' in path_B, do it
if [[ "$CLEAN_B" == "yes" || "$CLEAN_B" == "y" ]]; then
  echo "Running 'make clean' in $PATH_B..."
  cd "$PATH_B" || { echo "Failed to change directory to $PATH_B"; exit 1; }
  make clean || { echo "Make clean failed in $PATH_B"; exit 1; }
fi

# Now run the make target in the second directory (B)
echo "Running '$MAKE_TARGET' in $PATH_B..."
cd "$PATH_B" || { echo "Failed to change directory to $PATH_B"; exit 1; }
$MAKE_TARGET || { echo "Make failed in $PATH_B"; exit 1; }

echo "Makefiles ran successfully!"
