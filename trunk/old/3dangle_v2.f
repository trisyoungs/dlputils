!	** 3dangle v2 **
!	Program to calculate the angles (phi, theta) between a given atom (or the COM) of a molecule
!	and the defining X (phi) and Y (theta) axes of the molecule.
	IMPLICIT NONE
	EXTERNAL readframe,readheader,outinfo,calcgc,genaxis,getvector
	integer maxn,maxframes,maxbins,readframe,maxs,outinfo,patoms(2)
	integer readheader,maxm,grid,baselen
	PARAMETER(MAXN=8000,MAXFRAMES=1100,MAXBINS=400,MAXS=2,MAXM=200)
	PARAMETER(GRID=180)    ! Number of points in each direction, positive and negative, for grid
	character*80 header,infile,outfile,dlpoutfile,basename,resfile
	character*8 atmnam(MAXN)
	character*8 discard
	integer success,calctype,bin,X,numadded,comtype,n1,n2,n3,ngca(MAXS),gcal(MAXS,50),angtype
	integer angspecies
	integer keytrj, imcon, natms, nstep, iatm, n, discardn, nframes,k, l, fileform, m, o, p
	real*8 tstep, delta, px,py,pz, summ, vx,vy,vz,vmag,pi
	real*8 cell(9), xpos(MAXN),ypos(MAXN),zpos(MAXN), charge(MAXN), mass(MAXN)
	real*8 xvel(MAXN),yvel(MAXN),zvel(MAXN),xfor(MAXN),yfor(MAXN),zfor(MAXN)
	real*8 boxvolume, GCX(MAXS,MAXM),GCY(MAXS,MAXM),GCZ(MAXS,MAXM)
	real*8 tx,ty,tz,discardr,c1x,c1y,c1z,c2x,c2y,c2z,avcat(50,3)
	real*8 axisx(MAXS,MAXM,3),axisy(MAXS,MAXM,3),axisz(MAXS,MAXM,3)   ! Axis vectors (1:3 = x,y,z)
	real*8 ox(MAXS,MAXM),oy(MAXS,MAXM),oz(MAXS,MAXM)	    ! Axis origins
	real*8 angle(0:GRID,0:GRID),hist(2,0:GRID),a,b,c,cutoff
	real*8 sitex,sitey,sitez,rotx(10),roty(10),rotz(10),phi,theta,radcon
	real*8 v1x,v1y,v1z,v2x,v2y,v2z
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),aa(MAXS,3)
	COMMON /frame/ atmnam,xpos,ypos,zpos,charge,mass,xvel,yvel,zvel,xfor,yfor,zfor,
     +  cell,header,tstep,keytrj,imcon,natms,nstep
	COMMON /species/ s_name,s_natoms,s_nmols,s_start,nspecies
	COMMON /gc/ GCX,GCY,GCZ,gcal,comtype
	COMMON /axes/ axisx,axisy,axisz,ox,oy,oz

	pi=3.141592654
	comtype=1	! Type of centre-of-x calculations... (1=mass, 2=geometry)

	! Initialise some variables....
	DO n=0,GRID
	  DO m=0,GRID
	    angle(n,m)=0.0
	  END DO
	  hist(1,n)=0.0
	  hist(2,n)=0.0
	END DO

	WRITE(0,*) "*** 3Dangle"

	! Atom lists with which to generate the species axes
	aa(1,1)=1	! First two atoms specify x axis
	aa(1,2)=3
	aa(1,3)=2	! Third generates y

	! Create a list of atoms for use in COM calcs....
	DO n=1,8
	  gcal(1,n)=n
	END DO
	gcal(2,1)=1
	ngca(1)=8
	ngca(2)=1

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

	! Get the target centre....
	DO n=1,nspecies
	  WRITE(*,*) "Species ",n,":"
	  WRITE(*,"(8(I3,1X,'(',A5,')',3x))") (INT(X),atmnam(s_start(n)+X-1),X=1,s_natoms(n))
	END DO
	WRITE(0,*) "Enter the species that contains the centre:"
	READ(*,*) angspecies
	WRITE(0,*) "Enter the atom number, or '-1' for the COM:"
	READ(*,*) angtype
	WRITE(0,*) "Enter the distance cutoff to use:"
	READ(*,*) cutoff

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
	  CALL calcgc(n,n,xpos,ypos,zpos,cell,ngca(n),mass)
	  CALL genaxis(n,aa,xpos,ypos,zpos,cell)
	END DO
	! Loop over all molecules of species 1
	p=s_start(1)
	DO m=1,s_nmols(1)
	  numadded=0
	  ! Loop over all molecules of specified target species...
	  o=s_start(angspecies)
	  DO n=1,s_nmols(angspecies)
	    ! Grab the position that we're interested in... (use px,py,pz)
	    IF (angtype.NE.-1) THEN
	      px=xpos(o+angtype-1)
	      py=ypos(o+angtype-1)
	      pz=zpos(o+angtype-1)
	    ELSE
	      px=ox(1,n)
	      py=oy(1,n)
	      pz=oz(1,n)
	    END IF 
	    ! Calculate the vector between the origin of the molecule 'm' and this point....
	    v1x=px-xpos(p+7)    ! Unique H atom
	    v1y=py-ypos(p+7)
	    v1z=pz-zpos(p+7)
	    ! Normalise
	    vmag=SQRT(v1x*v1x + v1y*v1y + v1z*v1z)
	    !WRITE(0,*) vmag
	    ! Filter the species by distance
	    IF (vmag.LT.cutoff) THEN
	      v1x=v1x/vmag
	      v1y=v1y/vmag
	      v1z=v1z/vmag
	      ! Get reference vector (C -> Unique H vector)
	      v2x=xpos(p+7) - xpos(p+1)
	      v2y=ypos(p+7) - ypos(p+1)
	      v2z=zpos(p+7) - zpos(p+1)
	      vmag=SQRT(v2x*v2x + v2y*v2y + v2z*v2z)
	      v2x=v2x/vmag
	      v2y=v2y/vmag
	      v2z=v2z/vmag
	      ! Angle calculation
	      ! 1) Calculate angle between x-axis vector and v...
	      !phi=ACOS(vx*axisx(1,m,1) + vy*axisx(1,m,2) + vz*axisx(1,m,3))
	      phi=ACOS(v1x*v2x + v1y*v2y + v1z*v2z)
	      phi=(phi / PI) * 180.0
	      ! 2) Calculate angle between y-axis vector and v...
	      !theta=ACOS(vx*axisy(1,m,1) + vy*axisy(1,m,2) + vz*axisy(1,m,3))
	      theta=(theta / PI) * 180.0
	      ! Store the result.....
	      !WRITE(0,*) m,n,phi,theta,vmag
	      n1=NINT(phi)
	      n2=NINT(theta)
	      hist(1,n1)=hist(1,n1)+1
	      hist(2,n2)=hist(2,n2)+1
	      angle(n1,n2)=angle(n1,n2)+1
	    END IF
	    numadded=numadded+1
	    o=o+s_natoms(angspecies)
	  END DO   ! Loop over target molecules
	  p=p+s_natoms(1)
	END DO  ! Loop over cation molecules...
	! Next frame
	GOTO 101
120	WRITE(0,*) "Finished."
	! * Normalise the data *
	! **********************
	DO n1=0,GRID
	  DO n2=0,GRID
	    angle(n1,n2)=angle(n1,n2)/nframes  !/numadded
	  END DO
	  hist(1,n1)=hist(1,n1)  !/nframes/numadded
	  hist(2,n1)=hist(2,n1)  !/nframes/numadded
	END DO
	WRITE(0,*) "Check Value = ",angle(0,0)
	angle(0,0)=0.0
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
802	IF (baselen.EQ.-1) THEN
	  basename="3ddistresults."
	  baselen=14
	ELSE
	  basename=infile(1:baselen)
	ENDIF

	resfile=basename(1:baselen)//"angle"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	resfile=basename(1:baselen)//"histphi"
	OPEN(UNIT=10,file=resfile,FORM="FORMATTED")
	resfile=basename(1:baselen)//"histtheta"
	OPEN(UNIT=11,file=resfile,FORM="FORMATTED")
	DO n1=0,GRID
	  DO n2=0,GRID
	    WRITE(9,*) angle(n1,n2)
	  END DO
	  WRITE(10,*) n1,hist(1,n1)
	  WRITE(11,*) n1,hist(2,n1)
	END DO
	CLOSE(9)
	CLOSE(10)
	CLOSE(11)


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
	real*8 GCX(MAXS,MAXM),GCY(MAXS,MAXM),GCZ(MAXS,MAXM),discardr,massnorm
	real*8 x(MAXN),y(MAXN),z(MAXN),cell(9),tx,ty,tz,mass(MAXN),x0,y0,z0
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),gcal(MAXS,50)
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
	  GCX(s,n)=x(o+gcal(sp,1)-1)
	  GCY(s,n)=y(o+gcal(sp,1)-1)
	  GCZ(s,n)=z(o+gcal(sp,1)-1)
	  DO m=2,natms       ! Go through the atoms listed in gcal()
	    CALL PBC(x(o+gcal(sp,1)-1),y(o+gcal(sp,1)-1),z(o+gcal(sp,1)-1),GCX(s,n),GCY(s,n),GCZ(s,n),tx,ty,tz,cell)
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
	      CALL PBC(x(o+gcal(sp,m)-1),y(o+gcal(sp,m)-1),z(o+gcal(sp,m)-1),GCX(s,n),GCY(s,n),GCZ(s,n),tx,ty,tz,cell)
	      x0=x0 + (tx-GCX(s,n))*mass(o+gcal(sp,m)-1)
	      y0=y0 + (ty-GCY(s,n))*mass(o+gcal(sp,m)-1)
	      z0=z0 + (tz-GCZ(s,n))*mass(o+gcal(sp,m)-1)
	      massnorm=massnorm+mass(o+gcal(sp,m)-1)
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
	IF (s_natoms(sp).GT.2) THEN   ! Calculate local axes
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
	    ox(sp,n)=x(o+aa(sp,1)-1) + 0.5*vx
	    oy(sp,n)=y(o+aa(sp,1)-1) + 0.5*vy
	    oz(sp,n)=z(o+aa(sp,1)-1) + 0.5*vz
	    ! 3) Determine the Y axis of the molecule
	    !    -- Get the vector between aa atom 3 and midway along the x-axis vector
	    CALL PBC(x(o+aa(sp,3)-1),y(o+aa(sp,3)-1),z(o+aa(sp,3)-1),ox(sp,n),oy(sp,n),
     +    oz(sp,n),tx,ty,tz,cell)
	    CALL getvector(ox(sp,n),oy(sp,n),oz(sp,n),tx,ty,tz,vx,vy,vz)
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
	    !WRITE(88,*) "O: ",ox(sp,n),oy(sp,n),oz(sp,n)
	    !WRITE(88,*) "X: ",axisx(sp,n,1),axisx(sp,n,2),axisx(sp,n,3)
	    !WRITE(88,*) "Y: ",axisy(sp,n,1),axisy(sp,n,2),axisy(sp,n,3)
	    !WRITE(88,*) "Z: ",axisz(sp,n,1),axisz(sp,n,2),axisz(sp,n,3)
	  END DO
	ELSE
	  ! Otherwise, set the origin of the species to the coords of the first atom.
	  ! For, e.g., halide ions, this is a must!
	  o=s_start(sp)
	  DO n=1,s_nmols(sp)
	    ox(sp,n)=x(o)
	    oy(sp,n)=y(o)
	    oz(sp,n)=z(o)
	    o=o+s_natoms(sp)
	  END DO
	END IF
	END

	SUBROUTINE getvector(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	! xyz1,xyz2: coordinates of two atoms; xyz3: vector between them (result)
	real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3
	x3=x1-x2
	y3=y1-y2
	z3=z1-z2
	END
