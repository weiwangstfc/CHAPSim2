#!/bin/bash

# -----------------------------------------------------------------------------
# CHAPSim2 Build Script
# This script compiles the CHAPSim2 project while providing options for clean builds
# and debug mode configuration. It handles library dependencies and user interaction.
# -----------------------------------------------------------------------------

# Define relative paths to the directories containing the Makefiles
REL_PATH_LIB="./lib/2decomp-fft/build"
REL_PATH_BUILD="./build"
LIB_FILE="./lib/2decomp-fft/build/opt/lib/libdecomp2d.a"

# Resolve absolute paths from the relative paths
PATH_LIB=$(realpath "$REL_PATH_LIB" 2>/dev/null || readlink -f "$REL_PATH_LIB")
PATH_BUILD=$(realpath "$REL_PATH_BUILD" 2>/dev/null || readlink -f "$REL_PATH_BUILD")

# -----------------------------------------------------------------------------
# Step 1: Check and Build the Library
# -----------------------------------------------------------------------------
if [[ -f "$LIB_FILE" ]]; then
  echo "Library file '$LIB_FILE' exists. Skipping CMake in $PATH_LIB."
else
  echo "Library file not found. Running CMake in $PATH_LIB..."
  cd "$PATH_LIB" || { echo "Error: Failed to change directory to $PATH_LIB"; exit 1; }
  ./build_cmake.sh || { echo "Error: CMake build failed in $PATH_LIB"; exit 1; }
fi

# -----------------------------------------------------------------------------
# Step 2: Prompt for Build Options
# -----------------------------------------------------------------------------

# Prompt the user for 'make clean'
echo -n "Do you want to run 'make clean' in CHAPSim ($PATH_BUILD)? (yes/no): "
read CLEAN_BUILD

# Validate user input for 'make clean'
if [[ "$CLEAN_BUILD" != "yes" && "$CLEAN_BUILD" != "y" && "$CLEAN_BUILD" != "no" && "$CLEAN_BUILD" != "n" ]]; then
  echo "Invalid input. Please enter 'yes' or 'no'. Exiting."
  exit 1
fi

# Prompt the user for debug mode
echo -n "Do you want to run in debug mode? (yes/no): "
read DEBUG_MODE

# Validate debug mode input
if [[ "$DEBUG_MODE" != "yes" && "$DEBUG_MODE" != "y" && "$DEBUG_MODE" != "no" && "$DEBUG_MODE" != "n" ]]; then
  echo "Invalid input. Please enter 'yes' or 'no'. Exiting."
  exit 1
fi

# -----------------------------------------------------------------------------
# Step 3: Configure the Build Target
# -----------------------------------------------------------------------------

if [[ "$DEBUG_MODE" == "yes" || "$DEBUG_MODE" == "y" ]]; then
  echo "Debug mode enabled."
  if [[ "$CLEAN_BUILD" == "yes" || "$CLEAN_BUILD" == "y" ]]; then
    MAKE_TARGET="make all cfg=gnu"
  else
    MAKE_TARGET="make cfg=gnu"
  fi
else
  echo "Normal mode enabled."
  if [[ "$CLEAN_BUILD" == "yes" || "$CLEAN_BUILD" == "y" ]]; then
    MAKE_TARGET="make all"
  else
    MAKE_TARGET="make"
  fi
fi

# -----------------------------------------------------------------------------
# Step 4: Perform Clean Build if Requested
# -----------------------------------------------------------------------------

if [[ "$CLEAN_BUILD" == "yes" || "$CLEAN_BUILD" == "y" ]]; then
  echo "Running 'make clean' in $PATH_BUILD..."
  cd "$PATH_BUILD" || { echo "Error: Failed to change directory to $PATH_BUILD"; exit 1; }
  if ! make clean; then
    echo "Warning: 'make clean' failed in $PATH_BUILD. Continuing..."
  fi
fi

# -----------------------------------------------------------------------------
# Step 5: Run the Build
# -----------------------------------------------------------------------------

echo "Running '$MAKE_TARGET' in $PATH_BUILD..."
cd "$PATH_BUILD" || { echo "Error: Failed to change directory to $PATH_BUILD"; exit 1; }
if ! $MAKE_TARGET; then
  echo "Error: Build failed in $PATH_BUILD. Exiting."
  exit 1
fi

# -----------------------------------------------------------------------------
# Completion Message
# -----------------------------------------------------------------------------
echo "CHAPSim2 successfully compiled."