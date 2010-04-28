!	** rdf v2 **
!	RDF program to calculate pair-pair distribution functions between pairs of atom types
!	or a single atom type and the centre-of-mass of a second species.
!	Calculates the density at each frame of the simulation (for NPT simulations)
	EXTERNAL readframe,readheader,outinfo,mimdist,calcgc,calcpurerdf
	integer maxn,maxframes,maxbins,readframe,maxs,outinfo,gcal(50),patoms(2)
	integer readheader
	real*8 mimdist
	PARAMETER(MAXN=8000,MAXFRAMES=1100,MAXBINS=400,MAXS=20)
	character*80 header,infile,outfile,dlpoutfile
	character*8 atmnam(MAXN),speciesn(2)
	character*8 discard
	integer success,calctype,species(2),speciesa(2),bin,X,numadded,comtype
	integer keytrj, imcon, natms, nstep, iatm, n, discardn, nframes, l, fileform, m, o, p
	real*8 tstep, binwidth, dist
	real*8 cell(9), xpos(MAXN),ypos(MAXN),zpos(MAXN), charge(MAXN), mass(MAXN)
	real*8 xvel(MAXN),yvel(MAXN),zvel(MAXN),xfor(MAXN),yfor(MAXN),zfor(MAXN)
	real*8 boxvolume, GCX(2,300), GCY(2,300), GCZ(2,300)
	real*8 tx,ty,tz,discardr,c1x,c1y,c1z,c2x,c2y,c2z
	real*8 rdf(MAXBINS),purerdf(MAXBINS),nrdf(MAXBINS)
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),centre(2)
	COMMON /frame/ atmnam,xpos,ypos,zpos,charge,mass,xvel,yvel,zvel,xfor,yfor,zfor,
     +  cell,header,tstep,keytrj,imcon,natms,nstep
	COMMON /species/ s_name,s_natoms,s_nmols,s_start,nspecies
	COMMON /tcoords/ tx,ty,tz
	COMMON /gc/ GCX,GCY,GCZ,gcal,comtype
	COMMON /cpurerdf/ purerdf

	binwidth=0.1   ! In Angstroms
	DO n=1,MAXBINS
	  rdf(n)=0
	END DO

	! Type of centre-of-mass calculations...
	comtype=1	     ! 1 = centre of mass,     2 = geometric centre of mass
	WRITE(0,*) "*** rdf"
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
50	WRITE(0,*) "Name of the outfile to pull molecule data from? (0 for none)"
	READ(*,*) dlpoutfile
	IF (dlpoutfile(1:2).EQ."0 ") THEN  ! Get manual input.
	  WRITE(*,*) "Number of molecule types in config file?"
	  READ(*,*) nspecies
	  DO n=1,nspecies
	    WRITE(*,*) "Number of molecules of type: ",n,"?"
	    READ(*,*) s_nmols(n)
	    WRITE(*,*) "Number of atoms in molecule type: ",n,"?"
	    READ(*,*) s_natoms(n)
	  END DO
	ELSE     ! Read the data from an OUT file....
	  OPEN(UNIT=13,FILE=dlpoutfile,FORM="FORMATTED",IOSTAT=n)
	  IF (n.EQ.0) THEN   ! File opened successfully....
	    IF (outinfo(13,1).EQ.-1) GOTO 798
	    success=0
	  ELSE
	    success=-1
	  END IF
	END IF
	IF (success.EQ.-1) GOTO 50
	
	! Now, read in the atom names (i.e. the header) from the HISTORY file....
	IF (readheader(10,fileform).EQ.-1) GOTO 799

	! Get the information about the atom-pairs from the user now.....
	DO m=1,2
	  DO n=1,nspecies
	    WRITE(*,*) "Species ",n,":"
	    WRITE(*,"(8(I3,1X,'(',A5,')',3x))") (INT(X),atmnam(s_start(n)+X-1),X=1,s_natoms(n))
	  END DO
	  IF (m.EQ.1) WRITE(0,*) "Enter the type of the first centre:"
	  IF (m.EQ.2) WRITE(0,*) "Enter the type of the second centre:"
	  WRITE(0,*) "       1) Atom"
	  WRITE(0,*) "       2) Species (GC)       3) Species (Partial GC)"
55	  READ(*,*) centre(m)
	  IF (centre(m).EQ.1) THEN
	    WRITE(0,*) "Enter the species number which contains the atom:"
	    READ(*,*) species(m)
	    WRITE(0,*) "Enter the ID of the atom to use:"
	    READ(*,*) speciesa(m)
	  ELSE IF (centre(m).EQ.2) THEN
	    WRITE(0,*) "Enter the species number:"
	    READ(*,*) species(m)
	    patoms(m)=s_natoms(species(m))
	    DO n=1,patoms(m)
	      gcal(n)=n
	    END DO
	  ELSE IF (centre(m).EQ.3) THEN
	    WRITE(0,*) "Enter the species number:"
	    READ(*,*) species(m)
	    WRITE(0,*) "Select the subset of atoms in the species by:"
	    WRITE(0,*) "      1) Range"
	    WRITE(0,*) "      2) Explicit List"
	    READ(*,*) X
	    IF (X.EQ.1) THEN
	      WRITE(0,*) "First atom ID:"
	      READ(*,*) o
	      WRITE(0,*) "Second atom ID:"
	      READ(*,*) p
	      DO n=1,(p-o)+1
		gcal(n)=o+n-1
	      END DO
	      gcal(n+1)=-1   ! End the list
	      patoms(m)=(p-o)+1
	    ELSE IF (X.EQ.2) THEN
	      WRITE(0,*) "Enter atom IDs one at a time (-1 to end):"
	      n=0
51	      READ(*,*) o
	      IF (o.EQ.-1) THEN
		gcal(n+1)=-1
		patoms(m)=n
	      ELSE
		n=n+1
		gcal(n)=o
		GOTO 51
	      END IF
	    END IF
	  ELSE
	    GOTO 55
	  END IF
	  ! Now loop for second species...
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
	DO n=1,MAXBINS
	  rdf(n)=0
	END DO
	! Have coordinates of system.
	o=s_start(species(1))	 ! Set the start atom of the species
	! If necessary, calculate geometric centres now....
	IF (centre(1).GT.1) CALL calcgc(1,species(1),xpos,ypos,zpos,cell,patoms(1),mass)
	IF (centre(2).GT.1) CALL calcgc(2,species(2),xpos,ypos,zpos,cell,patoms(2),mass)
	! WRITE(0,*) (GCX(1,n),n=1,200)
	DO n=1,s_nmols(species(1))     ! Loop over all molecules of species 1...
	  IF (centre(1).EQ.1) THEN    ! Single atom
	    c1x=xpos(o+speciesa(1)-1)
	    c1y=ypos(o+speciesa(1)-1)
	    c1z=zpos(o+speciesa(1)-1)
	  ELSE     !  Grab the GC value grom the GCxyz arrays....
	    c1x=GCX(1,n)
	    c1y=GCY(1,n)
	    c1z=GCZ(1,n)
	  END IF
	  numadded=0
	  p=s_start(species(2))
	  ! Now loop over all molecules of second species....
	  DO m=1,s_nmols(species(2))
	    ! Set the coordinates for the second point....
	    IF (centre(2).EQ.1) THEN    ! Single atom
	      c2x=xpos(p+speciesa(2)-1)
	      c2y=ypos(p+speciesa(2)-1)
	      c2z=zpos(p+speciesa(2)-1)
	    ELSE     !  Grab the GC value grom the GCxyz arrays....
	      c2x=GCX(2,m)
	      c2y=GCY(2,m)
	      c2z=GCZ(2,m)
	    END IF
	    ! Get the shortest (MIM) distance between the atom pair...
	    CALL PBC(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz,cell)
	    dist=SQRT( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )
	    ! WRITE(0,*) c1x,c1y,c1z,c2x,c2y,c2z,tx,ty,tz,dist
	    ! 'Add' this distance...
	    bin=INT(dist*10)+1
	    rdf(bin)=rdf(bin)+1
	    numadded=numadded+1
	    p=p+s_natoms(species(2))
	  END DO
	  o=o+s_natoms(species(1))
	END DO
	! If the two species targets are *exactly* the same, we must decrease numadded by 1....
	! Calc and store the normalised rdf....
	CALL calcpurerdf(species(2),mass,binwidth,cell,0)
	DO n=1,MAXBINS
	  IF (centre(1).EQ.centre(2)) THEN
	    rdf(n)=rdf(n)/(numadded-1)
	  ELSE
	    rdf(n)=rdf(n)/numadded
	  END IF
	  nrdf(n)=nrdf(n)+rdf(n)/purerdf(n)
	END DO
	! Next frame
	GOTO 101
	! Work on the results now to get the proper RDF
120	WRITE(0,*) "Finished."
	! Also, we don't bother checkin if the two species are the same centres, so just
	! zero the first bin values....
	rdf(1)=0
	nrdf(1)=0
	  
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
	WRITE(0,*) "Name of the output file?"
	READ(*,*) outfile
	OPEN(UNIT=9,FILE=outfile,FORM="FORMATTED")
	WRITE(9,*) "Bin  N(Bin)    Pure     G(r)"
	DO n=1,MAXBINS
	  nrdf(n)=nrdf(n)/nframes
	END DO
	DO n=1,MAXBINS
	  WRITE(9,"(F4.1,3(',',F8.5))") n*binwidth,rdf(n),purerdf(n),nrdf(n)
	END DO
	CLOSE(9)
	WRITE(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	END

	REAL*8 FUNCTION mimdist(x1,y1,z1,x2,y2,z2,cell)
	! Performs MIM scheme. Parameters the two coordinate sets (NOT changed).
	REAL*8 x1,y1,z1,x2,y2,z2,tx,ty,tz,cell(9)
	!WRITE(0,*) "Original position:",x1,y1,z1
	!WRITE(0,*) "Refer. position:",x2,y2,z2
	CALL PBC(x1,y1,z1,x2,y2,z2,tx,ty,tz,cell)
	!WRITE(0,*) "MIMd position:",tx,ty,tz
	! Get the vector between the two points and return the distance
	mimdist=SQRT( (tx-x2)**2 + (ty-y2)**2 + (tz-z2)**2 )
	!WRITE(0,*) mimdist
	END

	SUBROUTINE PBC(x1,y1,z1,x2,y2,z2,x3,y3,z3,cell)
	REAL*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,cell(9)
	! Performs minimum image convention.....
	x3=x1
	y3=y1
	z3=z1
	IF (x1-x2.GT.cell(1)/2.0) x3=x3-cell(1)
	IF (x1-x2.LT.-cell(1)/2.0) x3=x3+cell(1)
	IF (y1-y2.GT.cell(5)/2.0) y3=y3-cell(5)
	IF (y1-y2.LT.-cell(5)/2.0) y3=y3+cell(5)
	IF (z1-z2.GT.cell(9)/2.0) z3=z3-cell(9)
	IF (z1-z2.LT.-cell(9)/2.0) z3=z3+cell(9)
	! New PBC
	x3=x1 - cell(1)*NINT((x1-x2)/cell(1))
	y3=y1 - cell(5)*NINT((y1-y2)/cell(5))
	z3=z1 - cell(9)*NINT((z1-z2)/cell(9))
	END

	SUBROUTINE calcgc(s,sp,x,y,z,cell,natms,mass)
	! Calculates the geometric/mass centres of all the molecules in species
	IMPLICIT NONE
	INTEGER MAXN,MAXS
	PARAMETER(MAXN=8000,MAXS=20)
	real*8 GCX(2,300),GCY(2,300),GCZ(2,300),discardr,massnorm
	real*8 x(MAXN),y(MAXN),z(MAXN),cell(9),tx,ty,tz,mass(MAXN),x0,y0,z0
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),gcal(50)
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
	  GCX(s,n)=x(o+gcal(1)-1)
	  GCY(s,n)=y(o+gcal(1)-1)
	  GCZ(s,n)=z(o+gcal(1)-1)
	  DO m=2,natms       ! Go through the atoms listed in gcal()
	    CALL PBC(GCX(s,n),GCY(s,n),GCZ(s,n),tx,ty,tz,tx,ty,tz,cell)
	    ! Add to the GCpos : multiply up first, then add, then divide for new GC
	    GCX(s,n)= ( (GCX(s,n)*(m-1)) + tx) / m
	    GCY(s,n)= ( (GCY(s,n)*(m-1)) + ty) / m
	    GCZ(s,n)= ( (GCZ(s,n)*(m-1)) + tz) / m
	  END DO
	  o=o+s_natoms(sp)
	END DO
	IF (comtype.EQ.1) THEN    ! Do centre-of-mass calculationa
	  o=s_start(sp)
	  DO n=1,s_nmols(sp)
	    ! We will use the geometric centre as the reference point....
	    x0=0.0
	    y0=0.0
	    z0=0.0
	    massnorm=0.0
	    DO m=1,natms
	      CALL PBC(x(o+gcal(m)-1),y(o+gcal(m)-1),z(o+gcal(m)-1),GCX(s,n),GCY(s,n),GCZ(s,n),tx,ty,tz,cell)
	      x0=x0 + (tx-GCX(s,n))*mass(o+gcal(m)-1)
	      y0=y0 + (ty-GCY(s,n))*mass(o+gcal(m)-1)
	      z0=z0 + (tz-GCZ(s,n))*mass(o+gcal(m)-1)
	      massnorm=massnorm+mass(o+gcal(m)-1)
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
	PARAMETER(MAXBINS=400,MAXN=8000,MAXS=20)
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
	
