#!/bin/bash

# -----------------------------------------------------------------------------
# CHAPSim2 Build Script
# This script compiles the CHAPSim2 project while providing options for clean builds
# and debug mode configuration. It handles library dependencies and user interaction.
# -----------------------------------------------------------------------------

# Define relative paths to the directories containing the Makefiles
REL_PATH_LIB="./lib/2decomp-fft/build"
REL_PATH_BUILD="./build"
REL_PATH_BIN="./bin"
CACHE_FILE="CMakeCache.txt"
LIB_FILE="$REL_PATH_LIB/opt/lib/libdecomp2d.a"

# Resolve absolute paths from the relative paths
PATH_LIB=$(realpath "$REL_PATH_LIB" 2>/dev/null || readlink -f "$REL_PATH_LIB")
PATH_BUILD=$(realpath "$REL_PATH_BUILD" 2>/dev/null || readlink -f "$REL_PATH_BUILD")
PATH_BIN=$(realpath "$REL_PATH_BIN" 2>/dev/null || readlink -f "$REL_PATH_BIN")

# Ensure the bin folder exists
if [[ ! -d "$PATH_BIN" ]]; then
  echo "Bin directory not found. Creating $PATH_BIN..."
  mkdir -p "$PATH_BIN" || { echo "Error: Failed to create bin directory"; exit 1; }
fi

# -----------------------------------------------------------------------------
# Step 1: Check and Build the Library
# -----------------------------------------------------------------------------
if [[ -f "$LIB_FILE" ]]; then
  echo "Library file '$LIB_FILE' exists. Skipping CMake in $PATH_LIB."
else
  echo "Library file not found. Running CMake in $PATH_LIB..."
  cd "$PATH_LIB" || { echo "Error: Failed to change directory to $PATH_LIB"; exit 1; }
  [[ -f "$CACHE_FILE" ]] && rm "$CACHE_FILE" || { echo "Error: Failed to remove $CACHE_FILE"; exit 1; }
  ./build_cmake.sh || { echo "Error: CMake build failed in $PATH_LIB"; exit 1; }
fi

# -----------------------------------------------------------------------------
# Step 2: Prompt for Build Options
# -----------------------------------------------------------------------------
read -p "Do you want to run 'make clean' before compiling CHAPSim? (yes/no/only): " CLEAN_BUILD

if [[ "$CLEAN_BUILD" != "yes" && "$CLEAN_BUILD" != "y" && "$CLEAN_BUILD" != "no" && "$CLEAN_BUILD" != "n" && "$CLEAN_BUILD" != "only" ]]; then
  echo "Invalid input. Please enter 'yes', 'no', or 'only'. Exiting."
  exit 1
fi

if [[ "$CLEAN_BUILD" == "only" ]]; then
  echo "Running 'make clean' in $PATH_BUILD..."
  cd "$PATH_BUILD" || { echo "Error: Failed to change directory to $PATH_BUILD"; exit 1; }
  if ! make clean; then
    echo "Warning: 'make clean' failed in $PATH_BUILD."
  fi
  exit 0
fi

read -p "Do you want to run in debug mode? (yes/no): " DEBUG_MODE

if [[ "$DEBUG_MODE" != "yes" && "$DEBUG_MODE" != "y" && "$DEBUG_MODE" != "no" && "$DEBUG_MODE" != "n" ]]; then
  echo "Invalid input. Please enter 'yes' or 'no'. Exiting."
  exit 1
fi

# -----------------------------------------------------------------------------
# Step 3: Configure the Build Target
# -----------------------------------------------------------------------------
if [[ "$DEBUG_MODE" =~ ^(yes|y)$ ]]; then
  echo "Debug mode enabled."
  MAKE_TARGET="make cfg=gnu"
else
  echo "Normal mode enabled."
  MAKE_TARGET="make"
fi

if [[ "$CLEAN_BUILD" =~ ^(yes|y)$ ]]; then
  MAKE_TARGET="make clean && $MAKE_TARGET all"
elif [[ "$CLEAN_BUILD" == "only" ]]; then
  MAKE_TARGET="make clean"
fi

# -----------------------------------------------------------------------------
# Step 4: Run the Build
# -----------------------------------------------------------------------------
echo "Running '$MAKE_TARGET' in $PATH_BUILD..."
cd "$PATH_BUILD" || { echo "Error: Failed to change directory to $PATH_BUILD"; exit 1; }
if ! eval $MAKE_TARGET; then
  echo "Error: Build failed in $PATH_BUILD. Exiting."
  exit 1
fi

# -----------------------------------------------------------------------------
# Completion Message
# -----------------------------------------------------------------------------
echo "CHAPSim2 successfully compiled."