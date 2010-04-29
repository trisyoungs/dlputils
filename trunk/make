#!/bin/bash

# Build all dlputils, setting up for the host machines listed below

HOST=`hostname`
INSTALL=false
OPTFLAGS="-O3"
for arg in "$*"  # Doesn't work properly if "$*" isn't quoted.
do
  if [ "$arg" = "install" ]; then
    INSTALL=true
    echo "... Will install build"
  elif [ "$arg" = "slow" ]; then
    OPTFLAGS="-CB -O0"
    echo "... Compiling with no optimisation"
  fi
done

if [ "$HOST" = "utopia" ]
then
  echo "Building on utopia..."
  FC=ifort
  FFLAGS="-i-static -static-libgcc $OPTFLAGS"
  MPIFC=ifort
  MPIFLAGS="-I../mpich-1.2.7p1/include -L../mpich-1.2.7p1/lib -lmpich"
elif [ "$HOST" = "n16" ]
then
  echo "Building on xcmaster..."
  FC=ifort
  FFLAGS="-i-static -static-libgcc $OPTFLAGS"
  MPIFC=/opt/hpmpi/bin/mpif90
  MPIFLAGS=""
fi

# Make
make FC=$FC FFLAGS="$FFLAGS" MPIFC=$MPIFC MPIFLAGS="$MPIFLAGS"

# Install?
if [ "$INSTALL" = "true" ]
then
  make install
fi

#FFLAGS=-ffree-line-length-180
