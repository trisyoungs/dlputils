!	** 3ddist v1 **
!	Program to calculate the 3-dimensional distribution of different species about
!	one species' centres.
	IMPLICIT NONE
	EXTERNAL readframe,readheader,outinfo,calcgc,genaxis,getvector
	integer maxn,maxframes,maxbins,readframe,maxs,outinfo,patoms(2)
	integer readheader,maxm,grid,baselen
	PARAMETER(MAXN=8000,MAXFRAMES=1100,MAXBINS=400,MAXS=2,MAXM=200)
	PARAMETER(GRID=20)    ! Number of points in each direction, positive and negative, for grid
	character*80 header,infile,outfile,dlpoutfile,basename,resfile
	character*8 atmnam(MAXN)
	character*8 discard
	integer success,calctype,bin,X,numadded,comtype,n1,n2,n3,ngca(MAXS),gcal(MAXS,50)
	integer keytrj, imcon, natms, nstep, iatm, n, discardn, nframes,k, l, fileform, m, o, p
	real*8 tstep, delta, px,py,pz, summ
	real*8 cell(9), xpos(MAXN),ypos(MAXN),zpos(MAXN), charge(MAXN), mass(MAXN)
	real*8 xvel(MAXN),yvel(MAXN),zvel(MAXN),xfor(MAXN),yfor(MAXN),zfor(MAXN)
	real*8 boxvolume, GCX(MAXS,MAXM),GCY(MAXS,MAXM),GCZ(MAXS,MAXM)
	real*8 tx,ty,tz,discardr,c1x,c1y,c1z,c2x,c2y,c2z,avcat(50,3)
	real*8 axisx(MAXS,MAXM,3),axisy(MAXS,MAXM,3),axisz(MAXS,MAXM,3)   ! Axis vectors (1:3 = x,y,z)
	real*8 ox(MAXS,MAXM),oy(MAXS,MAXM),oz(MAXS,MAXM)	    ! Axis origins
	real*8 pdens(MAXS,MAXS,-GRID:GRID,-GRID:GRID,-GRID:GRID)  ! Density distributions / number added
	real*8 ication(-GRID:GRID,-GRID:GRID,-GRID:GRID)         ! Intramolecular distribution for cation
	integer pdensn(MAXS,MAXS),pdensa(MAXS,MAXS)
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),aa(MAXS,3)
	COMMON /frame/ atmnam,xpos,ypos,zpos,charge,mass,xvel,yvel,zvel,xfor,yfor,zfor,
     +  cell,header,tstep,keytrj,imcon,natms,nstep
	COMMON /species/ s_name,s_natoms,s_nmols,s_start,nspecies
	COMMON /gc/ GCX,GCY,GCZ,gcal,comtype
	COMMON /axes/ axisx,axisy,axisz,ox,oy,oz

	delta=0.5 	! Spacing between grid points (Angstroms)
	comtype=1	! Type of centre-of-x calculations... (1=mass, 2=geometry)

	! Initialise some variables....
	DO l=1,MAXS
	  DO k=1,MAXS
	    pdensa(l,k) = 0
	pdensn(l,k) = 0
	    DO n=-GRID,GRID
	      DO m=-GRID,GRID
		DO o=-GRID,GRID
		  pdens(l,k,n,m,o)=0
		  ication(n,m,o)=0
		END DO
	      END DO
	    END DO
	  END DO
	END DO
	DO n=1,50
	  DO m=1,3
	    avcat(n,m)=0
	  END DO
	END DO

	WRITE(0,*) "*** 3Ddens"

	! Atom lists with which to generate the species axes
	aa(1,1)=1	! First two atoms specify x axis
	aa(1,2)=3
	aa(1,3)=2	! Third generates y

	! Create a list of atoms for use in COM calcs....
	DO n=1,10
	  gcal(1,n)=n
	END DO
	gcal(2,1)=1
	ngca(1)=10
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
	DO l=1,1	! Should be 1,nspecies but don't use anion-X data.
	  DO m=1,s_nmols(l)     ! Loop over all molecules of species 1...
	    ! Blank the pdensn() array here....
	    DO n=1,nspecies
	      !pdensn(l,n)=0
	    END DO
	    DO k=1,nspecies
	      p=s_start(k)
	      ! Now loop over all molecules of second species....
	      DO n=1,s_nmols(k)
		! * Calculations begin here *
		! ***************************
		! Loops are:
		!   l = Over all species 
		!     m = Over all molecules of species of l
		!       k = Over all species (2)
		!	 n = Over all molecules of species k
		! 1) Calculate the distribution of species k about species l
		!CALL PBC(GCX(k,n),GCY(k,n),GCZ(k,n),ox(l,m),oy(l,m),
		CALL PBC(ox(k,n),oy(k,n),oz(k,n),ox(l,m),oy(l,m),
     +		  oz(l,m),tx,ty,tz,cell)
		tx=tx-ox(l,m)
		ty=ty-oy(l,m)
		tz=tz-oz(l,m)
		!WRITE(40,*) GCX(k,n),GCY(k,n),GCZ(k,n),ox(l,m),oy(l,m),oz(l,m)
		!WRITE(40,*) tx,ty,tz
		! Apply a transformation to rotate into the axis of the molecule l
		px=tx*axisx(l,m,1) + ty*axisx(l,m,2) + tz*axisx(l,m,3)
		py=tx*axisy(l,m,1) + ty*axisy(l,m,2) + tz*axisy(l,m,3)
		pz=tx*axisz(l,m,1) + ty*axisz(l,m,2) + tz*axisz(l,m,3)
		!WRITE(40,*) px,py,pz
		! Use n1,n2,n3 to store the integers for the pdens() array
		n1=NINT(px/delta)
		n2=NINT(py/delta)
		n3=NINT(pz/delta)
		! If any of the n's are over GRID, then only increase the counter
		IF (MAX(ABS(n1),(MAX(ABS(n2),ABS(n3)))).GT.GRID) THEN
		!write(0,*) "REFUSED DIST",l,m,k,n,SQRT(tx**2 + ty**2 + tz**2)
		!write(0,*) "P",px,py,pz
		!write(0,*) "PDIST",SQRT(px**2 + py**2 + pz**2)
		!write(0,*) tx,ty,tz
		!write(0,*) n1,n2,n3
		  pdensn(l,k)=pdensn(l,k)+1      ! No position, just count it....
		ELSE
		  ! Check to make sure the same molecule isn't being consider with itself  XXXXX
		  pdens(l,k,n1,n2,n3)=pdens(l,k,n1,n2,n3)+1
		  pdensn(l,k)=pdensn(l,k)+1
		  pdensa(l,k)=pdensa(l,k)+1
		END IF
		p=p+s_natoms(k)
		! * Calculations END *
		! ********************
	      END DO
	    END DO
	    o=o+s_natoms(l)
	  END DO
	END DO
	! Calculate average cation and cation intramolecular distribution.
	p=s_start(1)
	DO n=1,s_nmols(1)
	  DO m=1,s_natoms(1)
	    ! PBC the atom
	    CALL PBC(xpos(p+m-1),ypos(p+m-1),zpos(p+m-1),ox(1,n),oy(1,n),oz(1,n),tx,ty,tz,cell)
	    ! 'Zero' its position with respect to the axis centre...
	    tx=tx-ox(1,n)
	    ty=ty-oy(1,n)
	    tz=tz-oz(1,n)
	    ! Transform the coordinate into the local coordinate system...
	    px=tx*axisx(1,n,1) + ty*axisx(1,n,2) + tz*axisx(1,n,3)
	    py=tx*axisy(1,n,1) + ty*axisy(1,n,2) + tz*axisy(1,n,3)
	    pz=tx*axisz(1,n,1) + ty*axisz(1,n,2) + tz*axisz(1,n,3)
	    ! Accumulate the position....
	    avcat(m,1)=avcat(m,1)+px
	    avcat(m,2)=avcat(m,2)+py
	    avcat(m,3)=avcat(m,3)+pz
	    ! Intramolecular distribution.....
	    n1=NINT(px/0.15)
	    n2=NINT(py/0.15)
	    n3=NINT(pz/0.15)
	    ication(n1,n2,n3)=ication(n1,n2,n3)+1
	  END DO
	  p=p+s_natoms(1)
	END DO
	! Next frame
	GOTO 101
120	WRITE(0,*) "Finished."
	! * Normalise the data *
	! **********************
	DO l=1,1	! Should be 1,nspecies but don't use anion-X data.
	  DO k=1,nspecies
	write(0,*) "PDENSN",l,k,pdensn(l,k)
	pdensn(l,k) = pdensn(l,k) / nframes / s_nmols(k)
	write(0,*) "  (Reduced = ",pdensn(l,k),")"
	   DO n1=-GRID,GRID
	    DO n2=-GRID,GRID
	     DO n3=-GRID,GRID
	      pdens(l,k,n1,n2,n3)=pdens(l,k,n1,n2,n3)/nframes/pdensn(l,k)/(delta**3)
	      ! pdens(l,k,n1,n2,n3)=pdens(l,k,n1,n2,n3)/nframes/200/(delta**3)
	     END DO
	    END DO
	   END DO
	  END DO
	END DO

	DO l=1,1
	  DO k=1,MAXS
	    write(0,*) "Original pdensa,n",l,k,pdensa(l,k),pdensn(l,k)
	write(0,*) s_nmols,nframes
	    pdensa(l,k)=pdensa(l,k)/nframes/s_nmols(k)/s_nmols(l)
	  END DO
	END DO
	DO n=1,s_natoms(1)
	  DO m=1,3
	    avcat(n,m)=avcat(n,m)/nframes/s_nmols(1)
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
801	WRITE(0,*) "Averages taken over ",nframes, " frames and ",pdensn(1,1)," atoms/species."
	WRITE(0,*) "Numerical density averages (caught within GRID boundaries):"
	DO n=1,1
	  DO m=1,MAXS
	    WRITE(0,*) "Species ",n," - Species ",m,": ",real(pdensa(n,m))/real((GRID*2+1)**3)
	  END DO
	END DO
	write(0,*) "Grid = ",(GRID*2+1),"x",(GRID*2+1),"x",(GRID*2+1)
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

	DO l=1,1	! Should be 1,nspecies, but anion is species 2 and we don't need that data.
	  DO k=1,MAXS
	    resfile=basename(1:baselen)//"3ddist"//CHAR(48+l)//CHAR(48+k)
	    OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	    DO n=-GRID,GRID
	      DO m=-GRID,GRID
		DO o=-GRID,GRID
		  ! WRITE(9,*) n*delta,m*delta,o*delta,pdens(l,k,n,m,o)
		  WRITE(9,*) pdens(l,k,n,m,o)
		END DO
	      END DO
	    END DO
	    CLOSE(9)
	  END DO
	END DO
	resfile=basename(1:baselen)//"avcat"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	! Write out a pdb of the average cation in a format that the CPMS module can understand
	DO n=1,s_natoms(1)
900	FORMAT ('ATOM  ',I5,A4,'0','UNK','I',I4,'R',3x,3F8.3)
	  WRITE(9,"('Atom',4x,I3,2x,A1,15x,3(F7.4,2x))") n,atmnam(s_start(1)-1+n)(1:1),avcat(n,1),avcat(n,2),avcat(n,3)
	!ATOM      1  C   ACh     1    5.4012   64.3027 61.5836
	END DO
	CLOSE(9)
	resfile=basename(1:baselen)//"icat"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	DO n=-GRID,GRID
	  DO m=-GRID,GRID
	    DO o=-GRID,GRID
	      WRITE(9,*) n*delta,m*delta,o*delta,ication(n,m,o)
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
