# CMake project file for dlputils

cmake_minimum_required (VERSION 2.6)
project (dlputils)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the NORMA directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -ffree-line-length-none")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, so will just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Setup module build
set (MODULES "parse" "dlprw" "utility" "pdensrw" "integerlist")
foreach (p ${MODULES})
  add_library(${p} "${p}.f90")
endforeach (p)

# setup binary builds
set (PROGRAMS "dlp2xyzf" "catdlp" "dlp2dlp" "addheader" "editconfig" "filterdlp" "his2config" "his2xyzs" "dlp2dlpflipz" "dlp2dlpshift" "dlpreorder" "chunkdlp" "xyz2his" "cdf" "dist2" "rdf" "rdfsum" "rdfaa" "rdfaainter" "rdfss" "rdfdep" "zdens" "acfcat" "acfprep" "acfsum" "angle" "cagecor2" "cn" "counthb" "geom" "getcell" "intertorsion" "intratorsion" "intratorsion2" "lifehb" "moldist" "msd" "msd2" "vac" "vac2" "bident" "bident2" "bident2anal" "bident3" "bident4" "cluster" "clusteroh" "cryscomp" "dahist" "legendre" "zangle" "glucanal" "glucgeom" "gluchb" "gluchbeach" "glucsphere" "pdens" "pdensint" "pdensgauss" "pdensmirror" "pdenstrim" "avgconfig" "dlpsize" "genlb" "gengg" "pairs" "probedlp")
foreach (p ${PROGRAMS})
  add_executable (${p} "${p}.f90")
  target_link_libraries(${p} "dlprw" "utility" "parse" "pdensrw" "integerlist")
endforeach (p)

# install executables and scripts
install (TARGETS ${PROGRAMS} 
         RUNTIME DESTINATION "bin")
