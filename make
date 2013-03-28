#!/bin/bash

# Custom options - modify these for 'unknown' machine builds
DEFAULTFC=gfortran
DEFAULTFFLAGS=" -O2 -ffree-line-length-none "
#DEFAULTFFLAGS=" -O2 "
#DEFAULTFFLAGS=" -CB "
DEFAULTMPIFC=mpif90
DEFAULTMPIFLAGS=""
OPTFLAGS="-O3 -g "
#OPTFLAGS="-CB -O0 -g "

# Build all dlputils, setting up for the host machines listed below

HOST=`hostname`
INSTALL=false
for arg in "$*"  # Doesn't work properly if "$*" isn't quoted.
do
  if [ "$arg" = "install" ]; then
    INSTALL=true
    echo "... Will install after build"
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
elif [ "$HOST" = "gorgon" ]
then
  echo "Building on gorgon..."
  FC=ifort
  FFLAGS="-i-static -static-libgcc $OPTFLAGS"
  MPIFC=ifort
  MPIFLAGS="-I../mpich-1.2.7p1/include -L../mpich-1.2.7p1/lib -lmpich"
else
  echo "Building on unknown machine..."
  FC=$DEFAULTFC
  FFLAGS=$DEFAULTFFLAGS
  MPIFC=$DEFAULTMPIFC
  MPIFLAGS=$DEFAULTMPIFLAGS
fi

# Make
make FC=$FC FFLAGS="$FFLAGS" MPIFC=$MPIFC MPIFLAGS="$MPIFLAGS"

# Install?
if [ "$INSTALL" = "true" ]
then
  make install
fi

