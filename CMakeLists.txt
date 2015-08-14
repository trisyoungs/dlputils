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

# Setup function library build
add_library(external "parse.f90" "dlprw.f90" "utility.f90" "pdensrw.f90" "integerlist.f90" "reallist.f90")

# Setup binary builds
set (PROGRAMS 
  "acfcat" "acfprep" "acfsum" "addheader" "angle" "ardf" "avgconfig"
  "bident" "bident2" "bident2anal" "bident3" "bident4" 
  "cagecor2" "catdlp" "cdf" "chunkdlp" "cluster" "clusterab" "cn" "counthb" "cryscomp"
  "dahist" "dist2" "dlp2dlp" "dlp2dlpflipz" "dlp2dlpshift" "dlp2xyzf" "dlpreorder"  "dlpsize"
  "editconfig" 
  "filterdlp"
  "getcell" "genlb" "gengg" "geom" "glucanal" "glucgeom" "gluchb" "gluchbeach"
  "his2config" "his2xyzs"
  "intertorsion" "intratorsion" "intratorsion2" 
  "legendre" "lifehb"
  "moldist" "msd" "msd2"
  "pairs" "pdens" "pdensgauss" "pdensint" "pdensmirror" "pdenstrim" "probedlp"  
  "rdf" "rdfaa" "rdfaainter" "rdfdep" "rdfss" "rdfsum" "ringconf"
  "xyz2his"
  "zangle" "zdens"
  "vac" "vac2"
)
foreach (p ${PROGRAMS})
  add_executable (${p} "${p}.f90")
  target_link_libraries(${p} "external")
endforeach (p)

# install executables and scripts
install (TARGETS ${PROGRAMS} 
         RUNTIME DESTINATION "bin")
