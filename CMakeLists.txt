# CMake project file for dlputils

cmake_minimum_required (VERSION 2.6)
project (dlputils)
set(DESCRIPTION "dlputils - Programs for analysis of DL_POLY trajectories")
set(AUTHOR "Tristan Youngs")
set(VERSION_MAJOR "1")
set(VERSION_MINOR "5")
set(VERSION_PATCH "0")

enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

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
add_library(external "parse.f90" "dlprw.f90" "utility.f90" "pdensrw.f90" "integerlist.f90" "reallist.f90" "ptable.f90")

# Setup binary builds
set (PROGRAMS 
  "acfcat" "acfprep" "acfsum" "addheader" "angle" "ardf" "avgconfig"
  "bident2" "bident2anal" "bident3" "bident4" 
  "cagecor2" "catdlp" "cdf" "cdf2" "chunkdlp" "cluster" "clusterab" "cn" "codf" "counthb" "cryscomp" "cdf_mol_map"
  "dahist" "dist2" "dlp2config" "dlp2dlp" "dlp2dlpflipz" "dlp2dlpshift" "dlp2xyzf" "dlp2xyzs" "dlpreorder"  "dlpsize" "dlp2xyz"
  "editconfig" 
  "dlpfilter"
  "getcell" "genlb" "gengg" "geom" "glucanal" "glucgeom" "gluchb" "gluchbeach"
  "intraangle" "intradist" "intertorsion" "intratorsion" "intratorsion2"
  "legendre" "lifehb"
  "moldist" "msd" "msd2"
  "orientrdf"
  "pairs" "pdens" "pdensgauss" "pdensinfo" "pdensint" "pdensmirror" "pdenstrim" "probedlp" "pdens_2d"  
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
         RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
