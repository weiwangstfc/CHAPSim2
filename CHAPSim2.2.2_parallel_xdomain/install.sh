#!/bin/sh

if make "$@"; then
  echo BUILD SUCCESSFUL
else
  echo BUILD FAILED
fi
