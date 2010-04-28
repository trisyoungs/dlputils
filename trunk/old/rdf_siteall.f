!	** rdf v5 **
!	RDF program to calculate pair-pair distribution functions between pairs of atom types
!	or a single atom type and the gemetric / mass centre of a second species.
!	Calculates the density at each frame of the simulation (for NPT simulations)
!	Assumes that all site-site g(r)'s between the cation (species1) atoms and the GC of the
!	second species is required.
	IMPLICIT NONE
	EXTERNAL readframe,readheader,outinfo,calcgc,calcpurerdf
	integer maxn,maxframes,maxbins,readframe,maxs,outinfo,maxla
	integer readheader
	PARAMETER(MAXN=8000,MAXFRAMES=1100,MAXBINS=400,MAXS=2,MAXLA=20)
	character*80 header,infile,outfile,dlpoutfile,basename,resfile
	character*8 atmnam(MAXN),speciesn(2)
	character*8 discard
	integer gcal(2,50),patoms(2),speciestype(MAXS),s1,s2,baselen
	integer success,calctype,bin,X,numadded,comtype,speciesa(MAXS)
	integer keytrj, imcon, natms, nstep, iatm, n, discardn, nframes, l, fileform, m, o, p
	real*8 tstep, binwidth, dist
	real*8 cell(9), xpos(MAXN),ypos(MAXN),zpos(MAXN), charge(MAXN), mass(MAXN)
	real*8 xvel(MAXN),yvel(MAXN),zvel(MAXN),xfor(MAXN),yfor(MAXN),zfor(MAXN)
	real*8 boxvolume, GCX(2,300), GCY(2,300), GCZ(2,300)
	real*8 tx,ty,tz,discardr,c1x,c1y,c1z,c2x,c2y,c2z
	real*8 rdf(MAXLA,MAXBINS),purerdf(MAXBINS),nrdf(MAXLA,MAXBINS)
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),centre(2)
	COMMON /frame/ atmnam,xpos,ypos,zpos,charge,mass,xvel,yvel,zvel,xfor,yfor,zfor,
     +  cell,header,tstep,keytrj,imcon,natms,nstep
	COMMON /species/ s_name,s_natoms,s_nmols,s_start,nspecies
	COMMON /gc/ GCX,GCY,GCZ,gcal,comtype
	COMMON /cpurerdf/ purerdf

	binwidth=0.1   ! In Angstroms
	DO m=1,MAXLA
	DO n=1,MAXBINS
	  rdf(m,n)=0.0
	END DO
	END DO

	! Type of centre-of-mass calculations...
	comtype=1	     ! 1 = centre of mass,     2 = geometric centre of mass
	WRITE(0,*) "*** rdf - Calculates all site-site[1],GC[2] functions"
	IF (comtype.EQ.1) WRITE(0,*) "( using true centres-of-mass for calculations )"
	IF (comtype.EQ.2) WRITE(0,*) "( using geometric centres for calculations )"

	WRITE(*,*) "Name of HISTORY file?"
	READ(*,*) infile
	! Determine the FORM of the HISTORY file...
	l=INDEX(infile," ")-1
	IF (infile(l:l).EQ."f") THEN
	  fileform=1
	ELSE IF (infile(l:l).EQ."F") THEN
	  fileform=1
	ELSE IF (infile(l:l).EQ."u") THEN
	  fileform=0
	ELSE IF (infile(l:l).EQ."U") THEN
	  fileform=0
	ELSE
	  WRITE(*,*) "Unable to determine the form of the HISTORY file."
	  WRITE(*,*) "Please enter 0 for unformatted or 1 for formatted:"
	  READ(*,*) fileform
	END IF

	IF (fileform.EQ.0) THEN
	  OPEN(UNIT=10,FILE=infile,FORM="UNFORMATTED")
	ELSE
	  OPEN(UNIT=10,FILE=infile,FORM="FORMATTED")
	END IF

	! Get the molecule data from the user....
	IF (outinfo(1).EQ.-1) GOTO 798
	
	! Now, read in the atom names (i.e. the header) from the HISTORY file....
	IF (readheader(10,fileform).EQ.-1) GOTO 799

	! No need to request info - we know what we're doing....
	patoms(2)=s_natoms(2)
	DO n=1,s_natoms(2)
	  gcal(2,n)=n
	END DO

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
100	nframes=0
	WRITE(0,*) "Fileform ",fileform
101	success=readframe(10,fileform)
	IF (success.EQ.1) GOTO 120  ! End of file encountered....
	IF (success.EQ.-1) GOTO 799  ! File error....
	nframes=nframes+1

	CALL calcgc(2,2,xpos,ypos,zpos,cell,patoms(2),mass)

	DO s1=1,s_natoms(1)     ! Loop over all atoms of species 1....

	DO n=1,MAXBINS
	  rdf(s1,n)=0
	END DO
	! Have coordinates of system.
	o=s_start(1)	 ! Set the start atom of the species
	DO n=1,s_nmols(1)     ! Loop over all molecules of species 1...
	  c1x=xpos(o+s1-1)
	  c1y=ypos(o+s1-1)
	  c1z=zpos(o+s1-1)
	  numadded=0
	  ! Now loop over all molecules of second species....
	  DO m=1,s_nmols(2)
	    ! Set the coordinates for the second point....
	    ! Grab the GC value grom the GCxyz arrays....
	    c2x=GCX(2,m)
	    c2y=GCY(2,m)
	    c2z=GCZ(2,m)
	    ! Get the shortest (MIM) distance between the atom pair...
	    CALL PBC(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz,cell)
	    dist=SQRT( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )
	    ! 'Add' this distance...
	    bin=INT(dist*10)+1
	    rdf(s1,bin)=rdf(s1,bin)+1
	    numadded=numadded+1
	  END DO
	  o=o+s_natoms(1)
	END DO

	END DO   ! End main loop over all atoms of species1.

	WRITE(0,*) "numadded:=",numadded
	! Calc and store the normalised rdf....
	CALL calcpurerdf(2,mass,binwidth,cell,0)
	DO s1=1,s_natoms(1)
	  DO n=1,MAXBINS
	    if (s1.EQ.1) write(0,*) rdf(s1,n),rdf(s1,n)/numadded
	    rdf(s1,n)=rdf(s1,n)/numadded
	    nrdf(s1,n)=nrdf(s1,n)+rdf(s1,n)/purerdf(n)
	  END DO
	END DO   ! End accumulation loop over all atoms of species 1.

	stop
	! Next frame
	GOTO 101
	! Work on the results now to get the proper RDF
120	WRITE(0,*) "Finished."
	GOTO 801

700	WRITE(*,*) "INFILE and OUTFILE have the same name!!!"
	GOTO 999

798	WRITE(0,*) "Problem with OUTPUT file."
	GOTO 999
799	WRITE(0,*) "HISTORY file ended prematurely!"
	GOTO 801
800	WRITE(0,*) "End of unformatted HISTORY file found."
801	WRITE(0,*) "Averages taken over ",nframes, " frames and ",numadded," atoms/species."
	WRITE(0,*) ""

	! Ascertain length of basename....
	baselen=-1
	DO n=80,1,-1
	  IF (infile(n:n).EQ.".") THEN
	    baselen=n
	    GOTO 802
	  ENDIF
	END DO
802     IF (baselen.EQ.-1) THEN
	  basename="rdfresults."
	  baselen=11
	ELSE
	  basename=infile(1:baselen)
	ENDIF

	
	DO s1=1,s_natoms(1)
	  resfile=basename(1:baselen)//"ssrdf"//CHAR(48+(s1/10))//CHAR(48+MOD(s1,10))//atmnam(s_start(1)+s1-1)(1:2)
	  OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	  !WRITE(9,*) "Bin     Raw     G(r)"
	  ! Normalise the RDFs with respect to the number of frames.
	  DO n=1,MAXBINS
	    nrdf(s1,n)=nrdf(s1,n)/nframes
	  END DO
	  DO n=1,MAXBINS
	    WRITE(9,"(F6.3,3x,F12.8)") (n*binwidth)-binwidth/2.0,nrdf(s1,n)
	  END DO
	  CLOSE(9)
	END DO
	WRITE(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	END

	SUBROUTINE PBC(x1,y1,z1,x2,y2,z2,x3,y3,z3,cell)
	! Performs minimum image convention.
	! xyz1: point to consider; xyz2: reference point; xyz3: minimum image coords of xyz1 (result)
	REAL*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,cell(9)
	! Performs minimum image convention.....
	x3=x1 - cell(1)*NINT((x1-x2)/cell(1))
	y3=y1 - cell(5)*NINT((y1-y2)/cell(5))
	z3=z1 - cell(9)*NINT((z1-z2)/cell(9))
	END

	SUBROUTINE calcgc(s,sp,x,y,z,cell,natms,mass)
	! Calculates the geometric/mass centres of all the molecules in species
	! : element of GCxyz array to fill; sp: species no; natms: no. of atoms in gcal()
	IMPLICIT NONE
	INTEGER MAXN,MAXS
	PARAMETER(MAXN=8000,MAXS=2)
	real*8 GCX(2,300),GCY(2,300),GCZ(2,300),discardr,massnorm
	real*8 x(MAXN),y(MAXN),z(MAXN),cell(9),tx,ty,tz,mass(MAXN),x0,y0,z0
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),gcal(2,50)
	integer s,sp,n,m,natms,o,p,comtype
	COMMON /gc/ GCX,GCY,GCZ,gcal,comtype
	COMMON /species/ s_name,s_natoms,s_nmols,s_start,nspecies
	DO n=1,s_nmols(sp)
	  GCX(s,n)=0.0
	  GCY(s,n)=0.0
	  GCZ(s,n)=0.0
	END DO
	! Always have to calculate the geometric centre first....
	o=s_start(sp)
	DO n=1,s_nmols(sp)
	  ! Add the first atom....
	  GCX(s,n)=x(o+gcal(s,1)-1)
	  GCY(s,n)=y(o+gcal(s,1)-1)
	  GCZ(s,n)=z(o+gcal(s,1)-1)
	  DO m=2,natms       ! Go through the atoms listed in gcal()
	    CALL PBC(x(o+gcal(s,1)-1),y(o+gcal(s,1)-1),z(o+gcal(s,1)-1),GCX(s,n),GCY(s,n),GCZ(s,n),tx,ty,tz,cell)
	    ! Add to the GCpos : multiply up first, then add, then divide for new GC
	    GCX(s,n)= ( (GCX(s,n)*(m-1)) + tx) / m
	    GCY(s,n)= ( (GCY(s,n)*(m-1)) + ty) / m
	    GCZ(s,n)= ( (GCZ(s,n)*(m-1)) + tz) / m
	  END DO
	  o=o+s_natoms(sp)
	END DO
	IF (comtype.EQ.1.AND.s_natoms(sp).GT.1) THEN    ! Do centre-of-mass calculation
	  o=s_start(sp)
	  DO n=1,s_nmols(sp)
	    ! We will use [0,0,0] as the reference point....
	    x0=0.0
	    y0=0.0
	    z0=0.0
	    massnorm=0.0
	    DO m=1,natms
	      CALL PBC(x(o+gcal(s,m)-1),y(o+gcal(s,m)-1),z(o+gcal(s,m)-1),GCX(s,n),GCY(s,n),GCZ(s,n),tx,ty,tz,cell)
	      x0=x0 + (tx-GCX(s,n))*mass(o+gcal(s,m)-1)
	      y0=y0 + (ty-GCY(s,n))*mass(o+gcal(s,m)-1)
	      z0=z0 + (tz-GCZ(s,n))*mass(o+gcal(s,m)-1)
	      massnorm=massnorm+mass(o+gcal(s,m)-1)
	    END DO
	    ! Normalise...
	    GCX(s,n)=GCX(s,n) + x0/massnorm
	    GCY(s,n)=GCY(s,n) + y0/massnorm
	    GCZ(s,n)=GCZ(s,n) + z0/massnorm 
	    o=o+s_natoms(sp)
	  END DO
	END IF
	END

	SUBROUTINE calcpurerdf(s,mass,width,cell,INFO)
	! Calculate the pure rdf
	IMPLICIT NONE
	INTEGER MAXBINS,INFO,MAXN,MAXS
	PARAMETER(MAXBINS=400,MAXN=8000,MAXS=2)
	real*8 species2m,species2d,species2v,purerdf(MAXBINS),mass(MAXN)
	real*8 const,pi,width,boxvolume,cell(9),widthcm
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),n,s
	COMMON /cpurerdf/ purerdf
	COMMON /species/ s_name,s_natoms,s_nmols,s_start,nspecies
	pi=3.141592654
	widthcm=width*1E-8
	! Work out the mass and density of species2.
	species2m=0
	DO n=1,s_natoms(s)
	  species2m=species2m+mass(s_start(s)+n-1)
	END DO
	IF(INFO.EQ.1) WRITE(0,*) "Calculated mass of species(2) = ",species2m," g/mol."
	species2m=species2m*s_nmols(s)
	species2m=species2m/6.02213E23     ! Actual g of species(2) in cell
	boxvolume=(cell(1)*cell(5)*cell(9))*(1E-8**3)   ! Volume in cm3
	species2d=species2m/boxvolume
	IF(INFO.EQ.1) WRITE(0,*) "Effective density of species(2) = ",species2d," g/cm3."
	species2v=(species2m/s_nmols(s))  
	IF(INFO.EQ.1) WRITE(0,*) "Molecular volume of species(2) = ",species2v," A3."
	const=(4.0*pi*species2d)/3.0
	DO n=1,MAXBINS
	  purerdf(n)=(const*((n*widthcm)**3-((n-1)*widthcm)**3))/species2v
	END DO
	RETURN
	END
