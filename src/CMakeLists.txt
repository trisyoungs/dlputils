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
