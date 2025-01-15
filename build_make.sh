#!/bin/bash

# Define the relative paths to the directories containing the Makefiles
REL_PATH_A="./lib/2decomp-fft/build"
REL_PATH_B="./build"
LIB_FILE="./lib/2decomp-fft/build/opt/lib/libdecomp2d.a"

# Get the absolute paths from the relative paths
PATH_A=$(realpath "$REL_PATH_A" 2>/dev/null || readlink -f "$REL_PATH_A")
PATH_B=$(realpath "$REL_PATH_B" 2>/dev/null || readlink -f "$REL_PATH_B")

# Check if the library file exists
if [[ -f "$LIB_FILE" ]]; then
  echo "$LIB_FILE exists. Skipping 'cmake' in $PATH_A."
else
  # Change to the first directory (A) and run the Makefile
  echo "Running CMake for libs in $PATH_A..."
  cd "$PATH_A" || { echo "Failed to change directory to $PATH_A"; exit 1; }
  ./build_cmake.sh || { echo "CMake build failed in $PATH_A"; exit 1; }
fi

# Ask the user if they want to run 'make clean' in PATH_B
echo "Do you want to run 'make clean' in CHAPSim $PATH_B? (yes/no): "
read CLEAN_B

# Validate user input for 'make clean'
if [[ "$CLEAN_B" != "yes" && "$CLEAN_B" != "y" && "$CLEAN_B" != "no" && "$CLEAN_B" != "n" ]]; then
  echo "Invalid input for 'make clean'. Please enter 'yes' or 'no'. Exiting."
  exit 1
fi

# Ask the user for the make target
echo "Enter the make target (e.g., 'make' or 'make all' or 'make all cfg=gnu'): "
read MAKE_TARGET

# Validate the make target
if ! [[ "$MAKE_TARGET" =~ ^make([[:space:]]+[-a-zA-Z0-9_=.]+)*$ ]]; then
  echo "Invalid make target: $MAKE_TARGET. Exiting."
  exit 1
fi

# If the user wants to run 'make clean' in PATH_B, do it
if [[ "$CLEAN_B" == "yes" || "$CLEAN_B" == "y" ]]; then
  echo "Running 'make clean' in $PATH_B..."
  cd "$PATH_B" || { echo "Failed to change directory to $PATH_B"; exit 1; }
  if ! make clean; then
    echo "Warning: Make clean failed in $PATH_B, continuing..."
  fi
fi

# Now run the make target in the second directory (B)
echo "Running '$MAKE_TARGET' in $PATH_B..."
cd "$PATH_B" || { echo "Failed to change directory to $PATH_B"; exit 1; }
$MAKE_TARGET || { echo "Make failed in $PATH_B"; exit 1; }

echo "Makefiles ran successfully!"