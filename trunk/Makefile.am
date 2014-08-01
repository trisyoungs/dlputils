bin_PROGRAMS = dlp2xyzf catdlp dlp2dlp addheader editconfig filterdlp his2config his2xyzs dlp2dlpflipz dlp2dlpshift dlpreorder chunkdlp xyz2his
bin_PROGRAMS += cdf dist2 rdf rdfsum rdfaa rdfaainter rdfss zdens
bin_PROGRAMS += acfcat acfprep acfsum angle cagecor2 cn counthb dlpgeom getcell intertorsion intratorsion lifehb moldist msd msd2 vac vac2
bin_PROGRAMS += bident bident2 bident2anal bident3 bident4 cluster cryscomp dahist legendre zangle
bin_PROGRAMS += glucanal glucgeom gluchb gluchbeach glucsphere
bin_PROGRAMS += pdens pdensint pdensmirror pdenstrim # pdensgauss
bin_PROGRAMS += avgconfig dlpsize genlb gengg pairs probedlp
bin_PROGRAMS += acf dacf sq totalgr

# Define some module variables
PARSE=parse.f90
DLPRW=dlprw.f90 
UTILITY=utility.f90
FIELDREAD=fieldread.f90

# In order to compile the last few programs, MPI is required.
# Set the following two lines to reflect the locations of the includes and libraries on your system
MPIFC="-I/usr/lib64/mpi/gcc/openmpi/include"
MPILD="-L/usr/lib64/mpi/gcc/openmpi/lib64 -lmpi_f90"

# Converters
addheader_SOURCES = ${PARSE} ${DLPRW} addheader.f90
catdlp_SOURCES = ${PARSE} ${DLPRW} catdlp.f90
chunkdlp_SOURCES = ${PARSE} ${DLPRW} chunkdlp.f90
dlp2xyzf_SOURCES = ${PARSE} ${DLPRW} dlp2xyzf.f90
dlp2dlp_SOURCES = ${PARSE} ${DLPRW} dlp2dlp.f90
dlp2dlpflipz_SOURCES = ${PARSE} ${DLPRW} dlp2dlpflipz.f90
dlp2dlpshift_SOURCES = ${PARSE} ${DLPRW} dlp2dlpshift.f90
dlpreorder_SOURCES = ${PARSE} ${DLPRW} dlpreorder.f90
editconfig_SOURCES = ${PARSE} ${DLPRW} editconfig.f90
filterdlp_SOURCES = ${PARSE} ${DLPRW} filterdlp.f90
his2config_SOURCES = ${PARSE} ${DLPRW} his2config.f90
his2xyzs_SOURCES = ${PARSE} ${DLPRW} his2xyzs.f90
xyz2his_SOURCES = ${PARSE} ${DLPRW} xyz2his.f90

# RDF
cdf_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} cdf.f90
dist2_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} dist2.f90
rdf_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} rdf.f90
rdfsum_SOURCES = ${PARSE} ${DLPRW} rdfsum.f90
rdfaa_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} rdf_aa.f90
rdfaainter_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} rdf_aa_inter.f90
rdfss_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} rdf_ss.f90
zdens_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} zdens.f90

# SQ
sq_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} sq.f90
sq_FCFLAGS = ${MPIFC}
sq_LDFLAGS = ${MPILD}
totalgr_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} totalgr.f90

# Analyse 1
acf_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} acf.f90
acf_FCFLAGS = ${MPIFC}
acf_LDFLAGS = ${MPILD}
acfcat_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} acfcat.f90
acfprep_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} acfprep.f90
acfsum_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} acfsum.f90
angle_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} angle.f90
cagecor2_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} cagecor2.f90
cn_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} cn.f90
counthb_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} counthb.f90
dacf_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} dacf.f90
dacf_FCFLAGS = ${MPIFC}
dacf_LDFLAGS = ${MPILD}
dlpgeom_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} dlpgeom.f90
getcell_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} getcell.f90
intertorsion_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} intertorsion.f90
intratorsion_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} intratorsion.f90
lifehb_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} lifehb.f90
moldist_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} moldist.f90
msd_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} msd.f90
msd2_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} msd2.f90
vac_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} vac.f90
vac2_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} vac2.f90

# Analyse 2
bident_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} bident.f90
bident2_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} bident2.f90
bident2anal_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} bident2anal.f90
bident3_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} bident3.f90
bident4_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} bident4.f90
cluster_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} cluster.f90
cryscomp_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} cryscomp.f90
dahist_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} dahist.f90
legendre_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} legendre.f90
zangle_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} zangle.f90

# Glucose
glucanal_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} glucanal.f90
glucgeom_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} glucgeom.f90
gluchb_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} gluchb.f90
gluchbeach_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} gluchbeach.f90
glucsphere_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} glucsphere.f90

# Probability Densities
pdens_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} pdens.f90
pdensint_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} pdensint.f90
pdensmirror_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} pdensmirror.f90
pdenstrim_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} pdenstrim.f90
#pdensgauss_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} pdensgauss.f90

# Misc
avgconfig_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} avgconfig.f90
dlpsize_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} dlpsize.f90
genlb_SOURCES = ${PARSE} genlb.f90
gengg_SOURCES = ${PARSE} gengg.f90
pairs_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} pairs.f90
probedlp_SOURCES = ${PARSE} ${DLPRW} ${UTILITY} probedlp.f90

