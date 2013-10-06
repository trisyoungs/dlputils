#=======================================================================
# Makefile for DL_POLY and MISC utils
#=======================================================================

#FC=gfortran
#FC=ifort
#FFLAGS=-ffree-line-length-180
#FFLAGS= -i-static -static-libgcc -CB -O0
#FFLAGS= -i-static -static-libgcc -O3
SHELL=/bin/bash

MODULES = parse.o dlprw.o utility.o fieldread.o
CONVERTERS = dlp2xyzf catdlp dlp2dlp addheader editconfig filterdlp his2config his2xyzs dlp2dlpflipz dlp2dlpshift dlpreorder chunkdlp xyz2his
RDF = rdf_ss rdf rdfsum rdf_aa rdf_aa_inter zdens dist2 cdf
SQ = totalgr sq 
ANALYSE = moldist dlpgeom counthb cn cagecor2 msd msd2 lifehb getcell angle intertorsion intratorsion vac vac2 acfsum acfprep acfcat dacf acf
ANALYSE2 = legendre cryscomp cluster zangle bident bident2 bident2anal bident3 bident4 dahist
GLUCOSE = glucgeom glucanal gluchb gluchbeach glucsphere
PDENS = pdens pdenstrim pdensmirror pdensint
MISC = dlpsize avgconfig pairs probedlp
MISCP = genlb gengg
ALL = $(CONVERTERS) $(RDF) $(ANALYSE2) $(GLUCOSE) $(PDENS) $(MISC) $(MISCP) $(ANALYSE) $(SQ)

all: $(MODULES) $(ALL)

convert: $(MODULES) $(CONVERTERS)
rdfs: $(MODULES) $(RDF)
analyse: $(MODULES) $(ANALYSE)
glucose: $(MODULES) $(GLUCOSE)
misc: $(MISC)

sq : sq.f90
	$(MPIFC) $(FFLAGS) -o $@ $< dlprw.o utility.o parse.o $(MPIFLAGS)

acf : acf.f90
	$(MPIFC) $(FFLAGS) -o $@ $< dlprw.o utility.o parse.o $(MPIFLAGS)

dacf : dacf.f90
	$(MPIFC) $(FFLAGS) -o $@ $< dlprw.o utility.o parse.o $(MPIFLAGS)

%.o : %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

% : %.f90
	$(FC) $(FFLAGS) -o $@ $< $(MODULES)

clean:
	/bin/rm *.o *.mod fort.* $(ALL)

install:
	rsync $(ALL) ~/bin
