!	** 3ddist v1 **
!	Program to calculate the 3-dimensional distribution of different species about
!	one species centres.
	IMPLICIT NONE
	EXTERNAL readframe,readheader,outinfo,calcgc,genaxis,getvector
	integer maxn,maxframes,maxbins,readframe,maxs,outinfo,gcal(50),patoms(2)
	integer readheader,maxm,grid
	PARAMETER(MAXN=8000,MAXFRAMES=1100,MAXBINS=400,MAXS=2,MAXM=200)
	PARAMETER(GRID=30)    ! Number of points in each direction, positive and negative, for grid
	character*80 header,infile,outfile,dlpoutfile,basename
	character*8 atmnam(MAXN),speciesn(2)
	character*8 discard
	integer success,calctype,species(2),speciesa(2),bin,X,numadded,comtype,n1,n2,n3
	integer keytrj, imcon, natms, nstep, iatm, n, discardn, nframes,k, l, fileform, m, o, p
	real*8 tstep, delta, px,py,pz
	real*8 cell(9), xpos(MAXN),ypos(MAXN),zpos(MAXN), charge(MAXN), mass(MAXN)
	real*8 xvel(MAXN),yvel(MAXN),zvel(MAXN),xfor(MAXN),yfor(MAXN),zfor(MAXN)
	real*8 boxvolume, GCX(2,MAXM), GCY(2,MAXM), GCZ(2,MAXM)
	real*8 tx,ty,tz,discardr,c1x,c1y,c1z,c2x,c2y,c2z
	real*8 axisx(MAXS,MAXM,3),axisy(MAXS,MAXM,3),axisz(MAXS,MAXM,3)   ! Axis vectors (1:3 = x,y,z)
	real*8 ox(MAXS,MAXM),oy(MAXS,MAXM),oz(MAXS,MAXM)	    ! Axis origins
	real*8 pdens(MAXS,MAXS,GRID,GRID,GRID)  ! Density distributions / number added
	integer pdensn(MAXS,MAXS)
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),centre(2),aa(MAXS,3)
	COMMON /frame/ atmnam,xpos,ypos,zpos,charge,mass,xvel,yvel,zvel,xfor,yfor,zfor,
     +  cell,header,tstep,keytrj,imcon,natms,nstep
	COMMON /species/ s_name,s_natoms,s_nmols,s_start,nspecies
	COMMON /gc/ GCX,GCY,GCZ,gcal,comtype
	COMMON /axes/ axisx,axisy,axisz,ox,oy,oz

	delta=0.5	! Spacing between grid points (Angstroms)
	comtype=1	! Type of centre-of-x calculations... (1=mass, 2=geometry)

	! Initialise some variables....
	DO l=1,MAXS
	  DO k=1,MAXS
	    DO n=1,GRID
	      DO m=1,GRID
		DO o=1,GRID
		  pdens(l,k,n,m,o)=0
		END DO
	      END DO
	    END DO
	  END DO
	END DO

	WRITE(0,*) "*** 3Ddens"
	IF (comtype.EQ.1) WRITE(0,*) "( using true centres-of-mass for calculations )"
	IF (comtype.EQ.2) WRITE(0,*) "( using geometric centres for calculations )"

	! Atom lists with which to generate the species axes
	aa(1,1)=1	! First two atoms specify x axis
	aa(1,2)=3
	aa(1,3)=2	! Third generates y

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

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
	WRITE(0,*) "Fileform ",fileform
101	success=readframe(10,fileform)
	IF (success.EQ.1) GOTO 120  ! End of file encountered....
	IF (success.EQ.-1) GOTO 799  ! File error....
	nframes=nframes+1
	WRITE(0,*) nframes
	! Calculate all geometric centres....
	DO n=1,nspecies
	  CALL calcgc(n,n,xpos,ypos,zpos,cell,patoms(n),mass)
	  IF (s_natoms(n).GT.2) CALL genaxis(n,aa,xpos,ypos,zpos,cell)
	END DO
	DO l=1,nspecies
	  DO m=1,s_nmols(species(l))     ! Loop over all molecules of species 1...
	    IF (centre(1).EQ.1) THEN    ! Single atom
	      c1x=xpos(o+speciesa(l)-1)
	      c1y=ypos(o+speciesa(l)-1)
	      c1z=zpos(o+speciesa(l)-1)
	    ELSE     !  Grab the GC value grom the GCxyz arrays....
	      c1x=GCX(l,m)
	      c1y=GCY(l,m)
	      c1z=GCZ(l,m)
	    END IF
	    numadded=0
	    DO k=1,nspecies
	      p=s_start(species(k))
	      ! Now loop over all molecules of second species....
	      DO n=1,s_nmols(species(k))
		! Set the coordinates for the second point....
		IF (centre(k).EQ.1) THEN    ! Single atom
		  c2x=xpos(p+speciesa(k)-1)
		  c2y=ypos(p+speciesa(k)-1)
		  c2z=zpos(p+speciesa(k)-1)
		ELSE     !  Grab the GC value grom the GCxyz arrays....
		  c2x=GCX(k,n)
		  c2y=GCY(k,n)
		  c2z=GCZ(k,n)
		END IF
		! * Calculations begin here *
		! ***************************
		! Loops are:
		!   l = Over all species 
		!     m = Over all molecules of species of l
		!       k = Over all species (2)
		!	 n = Over all molecules of species k
		! 1) Calculate the distribution of species k about species l
		CALL PBC(GCX(k,n),GCY(k,n),GCZ(k,n),GCX(l,m),GCY(l,m),
     +		  GCZ(l,m),tx,ty,tz,cell)
		tx=tx-GCX(l,m)
		ty=ty-GCY(l,m)
		tz=tz-GCZ(l,m)
		! Apply a transformation to rotate into the axis of the molecule l
		px=tx*axisx(l,m,1) + ty*axisx(l,m,2) + tz*axisx(l,m,3)
		py=tx*axisy(l,m,1) + ty*axisy(l,m,2) + tz*axisy(l,m,3)
		pz=tx*axisz(l,m,1) + ty*axisz(l,m,2) + tz*axisz(l,m,3)
		! Use tx,ty,tz to store the integers for the pdens() array
		n1=NINT(px/delta)
		n2=NINT(py/delta)
		n3=NINT(pz/delta)
		! If any of the n's are over GRID/2, then only increase the counter
		IF (MAX(n1,(MAX(n2,n3))).GT.GRID/2) THEN
		  pdensn(l,k)=pdensn(l,k)+1      ! No position, just count it....
		ELSE
		  ! Check to make sure the same molecule isn't being consider with itself
		  pdens(l,k,n1,n2,n3)=pdens(l,k,n1,n2,n3)+1
		  pdensn(l,k)=pdensn(l,k)+1
		END IF
		p=p+s_natoms(species(k))
		! * Calculations END *
		! ********************
	      END DO
	    END DO
	    o=o+s_natoms(species(l))
	  END DO
	END DO
	! Next frame
	GOTO 101
120	WRITE(0,*) "Finished."
	! * Normalise the data *
	! **********************
	DO l=1,nspecies
	  DO k=1,nspecies
	   DO n1=1,GRID
	    DO n2=1,GRID
	     DO n3=1,GRID
	      pdens(l,k,n1,n2,n3)=pdens(l,k,n1,n2,n3)/nframes/pdensn(l,k)/(delta**3)
	     END DO
	    END DO
	   END DO
	  END DO
	END DO
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
	basename=infile(1:n)
	! Check
	IF (baselen.EQ.-1) basename="3ddistresults."

	! Break up output into individual species, and limit to a range of 15:15
	DO l=1,MAXS
	  DO k=1,MAXS
	    resfile=basename//CHAR(48+l)//CHAR(48+k)
	    OPEN(UNIT=9,file=resfile,FORM=FORMATTED)
	    DO n=1,GRID
	      DO m=1,GRID
		DO o=1,GRID
		  WRITE(9,*) n*delta,m*delta,o*delta,pdens(l,k,n,m,o)
		END DO
	      END DO
	    END DO
	  END DO
	END DO

	CLOSE(9)
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
	! s: element of GCxyz array to fill; sp: species no; natms: no. of atoms in gcal()
	IMPLICIT NONE
	INTEGER MAXN,MAXS,MAXM
	PARAMETER(MAXN=8000,MAXS=2,MAXM=200)
	real*8 GCX(2,MAXM),GCY(2,MAXM),GCZ(2,MAXM),discardr,massnorm
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
	    CALL PBC(x(o+gcal(1)-1),y(o+gcal(1)-1),z(o+gcal(1)-1),GCX(s,n),GCY(s,n),GCZ(s,n),tx,ty,tz,cell)
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

	SUBROUTINE genaxis(sp,aa,x,y,z,cell)
	! Subroutine to calculate the common axis of the species
	INTEGER MAXN,MAXS,MAXM
	PARAMETER(MAXN=8000,MAXS=2,MAXM=200)
	real*8 vx,vy,vz,tx,ty,tz,ox(MAXS,MAXM),oy(MAXS,MAXM),oz(MAXS,MAXM)
	real*8 axisx(MAXS,MAXM,3),axisy(MAXS,MAXM,3),axisz(MAXS,MAXM,3),vmag
	real*8 x(MAXN),y(MAXN),z(MAXN),cell(9)
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),sp,aa(MAXS,3)
	COMMON /species/ s_name,s_natoms,s_nmols,s_start,nspecies
	COMMON /axes/ axisx,axisy,axisz,ox,oy,oz
	o=s_start(sp)
	DO n=1,s_nmols(sp)
	  ! 1) Determine the X axis of the molecule
	  !    -- Get the vector between atoms 1 and 2...
	  CALL PBC(x(o+aa(sp,2)-1),y(o+aa(sp,2)-1),z(o+aa(sp,2)-1),
     +    x(o+aa(sp,1)-1),y(o+aa(sp,1)-1),z(o+aa(sp,1)-1),tx,ty,tz,cell)
	  CALL getvector(tx,ty,tz,x(o+aa(sp,1)-1),y(o+aa(sp,1)-1),z(o+aa(sp,1)-1),vx,vy,vz)
	  vmag=SQRT(vx**2 + vy**2 + vz**2)     ! Magnitude of vector
	  axisx(sp,n,1)=vx/vmag	! Store the normalised components of the x-axis vector
	  axisx(sp,n,2)=vy/vmag
	  axisx(sp,n,3)=vz/vmag
	  ! 2) Set the origin of the axis system - halfway along the x-axis vector
	  ox(sp,n)=x(o+aa(sp,1)-1) + 0.5*vx*vmag
	  oy(sp,n)=y(o+aa(sp,1)-1) + 0.5*vy*vmag
	  oz(sp,n)=z(o+aa(sp,1)-1) + 0.5*vz*vmag
	  ! 3) Determine the Y axis of the molecule
	  !    -- Get the vector between atom 3 and midway along the x-axis vector
	  CALL PBC(x(o+aa(sp,3)-1),y(o+aa(sp,3)-1),z(o+aa(sp,3)-1),x(o+aa(sp,1)-1)+0.5*vx,
     +    y(o+aa(sp,1)-1)+0.5*vy,z(o+aa(sp,1)-1)+0.5*vz,tx,ty,tz,cell)
	  CALL getvector(tx,ty,tz,ox(sp,n),oy(sp,n),oz(sp,n),vx,vy,vz)
	  vmag=SQRT(vx**2 + vy**2 + vz**2)     ! Magnitude of vector
	  axisy(sp,n,1)=vx/vmag   ! Store the normalised components of the y-axis vector
	  axisy(sp,n,2)=vy/vmag
	  axisy(sp,n,3)=vz/vmag
	  ! 4) Calculate the z axis from the x and y axes....
	  axisz(sp,n,1)=axisx(sp,n,2)*axisy(sp,n,3) - axisx(sp,n,3)*axisy(sp,n,2)
	  axisz(sp,n,2)=axisx(sp,n,3)*axisy(sp,n,1) - axisx(sp,n,1)*axisy(sp,n,3)
	  axisz(sp,n,3)=axisx(sp,n,1)*axisy(sp,n,2) - axisx(sp,n,2)*axisy(sp,n,1)
	  ! Increase the atom count for the next molecule....
	  o=o+s_natoms(sp)
	END DO
	END

	SUBROUTINE getvector(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	! xyz1,xyz2: coordinates of two atoms; xyz3: vector between them (result)
	real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3
	x3=x2-x1
	y3=y2-y1
	z3=z2-z1
	END
