#=======================================================================
# Makefile for DL_POLY and MISC utils
#=======================================================================

#
#FC=gfortran
FC=ifort
#FFLAGS=-ffree-line-length-180
#FFLAGS= -i-static -static-libgcc -CB -O0
FFLAGS= -i-static -static-libgcc -O3
SHELL=/bin/bash

MODULES = parse.o dlprw.o utility.o fieldread.o
CONVERTERS = dlp2fm catdlp config2bgf dlp2dlp addheader editconfig everydlp filterdlp his2config his2xyzs thindlp dlp2dlpflipz dlp2dlpshift dlpreorder
RDF = rdf_ss rdf sk rdfsum rdf_aa rdf_aa_inter zdens dist2 raw2sq
SQ = totalgr sk sq sqtest
ANALYSE = moldist dlpgeom counthb counthb_il cn cagecor msd msd2 lifehb getcell angle intertorsion intratorsion vac vac2
ANALYSE2 = axishist axishist2 legendre cryscomp cluster zangle bident bident2 bident2anal dahist #planedist
GLUCOSE = glucgeom glucanal gluchb gluchbeach glucsphere
PDENS = pdens pdenstrim pdensmirror pdensint rpairs surfacify
MISC = prepfm dlpsize point avgconfig pairs probedlp velconfig
MISCP = genlj 
ALL = $(MODULES) $(CONVERTERS) $(RDF) $(ANALYSE2) $(ANALYSE) $(GLUCOSE) $(PDENS) $(MISC) $(MISCP) $(SQ)

all: $(ALL)

convert: $(MODULES) $(CONVERTERS)
rdfs: $(MODULES) $(RDF)
analyse: $(MODULES) $(ANALYSE)
glucose: $(MODULES) $(GLUCOSE)
misc: $(MISC) parse.o

sq : sq.f90
	$(FC) $(FFLAGS) -o $@ $< dlprw.o utility.o parse.o -L/home/tris/src/mpich-1.2.7p1/lib -lmpich

sqtest : sqtest.f90
	$(FC) $(FFLAGS) -o $@ $< dlprw.o utility.o parse.o -L/home/tris/src/mpich-1.2.7p1/lib -lmpich

%.o : %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

% : %.f90
	$(FC) $(FFLAGS) -o $@ $< $(MODULES)

clean:
	/bin/rm *.o *.mod fort.* $(ALL)

install:
	rsync $(ALL) ~/bin
