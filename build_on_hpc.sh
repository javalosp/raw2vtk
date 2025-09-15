#!/usr/bin/env bash

# Clean previous build
make clean

# Create release directory if it doesn't exist
mkdir -p release

# This is just to make sure that we are using the appropriate version of modules
module purge

## Load system modules
module load vtk/5.8.0 
#module load mpi/intel-2019
module load mpi/intel-2019.8.254 # This is prefered over the default 2019
# CX3 versions
#module load tools/prod
#module load iimpi/2021a
# Using CX3 version of mpi compilers makes boost unhappy
# Therefore we need to load intel suite
#module load intel-suite/2019.4
module load boost/1.72.0 # This also loads intel-suite/2019.4 (boost requirement)

make uint8 VERBOSE=1
