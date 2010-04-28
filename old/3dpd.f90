!	** 3dpd.f90 **
!	  (based on old 3ddist_v2)
!	Program to calculate the 3-dimensional distribution of different species about
!	one species' centres.

	module pdens_data
	  real*8, allocatable :: gcx(:,:), gcy(:,:), gcz(:,:), pdens(:,:,:,:)
	  real*8, allocatable :: axisx(:,:,:), axisy(:,:,:), axisz(:,:,:)   ! Axis point coords
	  real*8, allocatable :: axisox(:,:), axisoy(:,:), axisoz(:,:)	    ! Axis origin
	  integer, allocatable :: pdensn(:), pdensa(:)
	end module pdens_data

	program pdens
	use dlprw
	implicit none
	integer :: maxframes,maxbins,patoms(2)
	integer :: baselen
	PARAMETER(MAXFRAMES=1100,MAXBINS=400)
	character*80 :: hisfile,outfile,basename,resfile
	character*8 :: discard
	integer :: success,calctype,bin,X,numadded,comtype,n1,n2,n3
	integer :: iatm, n, discardn, nframes,k, l, m, o, p
	real*8  delta, px,py,pz, summ
	real*8 boxvolume
	real*8 tx,ty,tz,discardr,c1x,c1y,c1z,c2x,c2y,c2z,avcat(50,3)
	!real*8 pdens(MAXS,MAXS,-GRID:GRID,-GRID:GRID,-GRID:GRID)  ! Density distributions / number added
	real*8 ication(-GRID:GRID,-GRID:GRID,-GRID:GRID)         ! Intramolecular distribution for cation
	integer :: aa(MAXS,3),centresp,sp1,m1,sp2,m2

	delta=0.5 	! Spacing between grid points (Angstroms)
	comtype=1	! Type of centre-of-x calculations... (1=mass, 2=geometry)
	grid = 20	! Number of 'boxes' in each -x,x,-y,y,-z,z direction

	write(0,*) "*** 3Dpd"

	! Atom lists with which to generate the species axes
	aa(1,1)=1	! First two atoms specify x axis
	aa(1,2)=3
	aa(1,3)=2	! Third generates y

	write(*,*) "Name of HISTORY file?"
	READ(*,*) hisfile
	call openhis(hisfile,10)
	! Get the molecule data from the user....
	write(0,*) "OUT File?"
	read(*,*) outfile
	IF (outinfo(outfile,1).EQ.-1) GOTO 798

	call alloc_xyz
	call alloc_vars
	
	! Now, read in the atom names (i.e. the header) from the HISTORY file....
	IF (readheader().EQ.-1) GOTO 799

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	IF (success.EQ.1) GOTO 120  ! End of file encountered....
	IF (success.EQ.-1) GOTO 799  ! File error....
	nframes=nframes+1
	write(0,*) nframes
	! Calculate all geometric centres....
	call calcgc
	! Generate all molecular axes
	call genaxis(n,aa,xpos,ypos,zpos,cell)


	sp1 = centremol
	do m1=1,s_nmols(sp1)     ! Loop over all molecules of species 1...
	  do sp2=1,nspecies
	    ! Now loop over all molecules of second species....
	    do n=1,s_nmols(sp2)
	      ! * Calculations begin here *
	      ! ***************************
	      ! 1) Calculate the distribution of species k about species l
	      call pbc(ox(sp2,m2),oy(sp2,m2),oz(sp2,m2),ox(sp1,m1),oy(sp1,m1), &
		& oz(sp1,m1),tx,ty,tz)
	      tx=tx-ox(sp1,m1)
	      ty=ty-oy(sp1,m1)
	      tz=tz-oz(sp1,m1)
	      !write(40,*) gcx(k,n),gcy(k,n),gcz(k,n),ox(l,m),oy(l,m),oz(l,m)
	      !write(40,*) tx,ty,tz
	      ! Apply a transformation to rotate into the axis of the molecule l
	      px=tx*axisx(sp1,m1,1) + ty*axisx(sp1,m1,2) + tz*axisx(sp1,m1,3)
	      py=tx*axisy(sp1,m1,1) + ty*axisy(sp1,m1,2) + tz*axisy(sp1,m1,3)
	      pz=tx*axisz(sp1,m1,1) + ty*axisz(sp1,m1,2) + tz*axisz(sp1,m1,3)
	      !write(40,*) px,py,pz
	      ! Use n1,n2,n3 to store the integers for the pdens() array
	      n1=NINT(px/delta)
	      n2=NINT(py/delta)
	      n3=NINT(pz/delta)
	      ! If any of the n's are over GRID, then only increase the counter
	      if (MAX(ABS(n1),(MAX(ABS(n2),ABS(n3)))).GT.GRID) THEN
	        write(0,*) "REFUSED DIST",l,m,k,n,SQRT(tx**2 + ty**2 + tz**2)
	        write(0,*) "P",px,py,pz
	        write(0,*) "PDIST",SQRT(px**2 + py**2 + pz**2)
	        write(0,*) tx,ty,tz
	        write(0,*) n1,n2,n3
	        pdensn(sp2)=pdensn(sp2)+1      ! No position, just count it....
	      else
	        ! Check to make sure the same molecule isn't being consider with itself  XXXXX
	        pdens(sp2,n1,n2,n3)=pdens(sp2,n1,n2,n3)+1
	        pdensn(sp2)=pdensn(sp2)+1
	        pdensa(sp2)=pdensa(sp2)+1
	      end if
	      p=p+s_natoms(sp2)
	      ! * Calculations end *
	      ! ********************
	    end do
	  end do
	end do
	! Calculate average cation and cation intramolecular distribution.

	p=s_start(sp1)
	do m1=1,s_nmols(sp1)
	  do n=1,s_natoms(sp1)
	    ! PBC the atom
	    CALL PBC(xpos(p+n-1),ypos(p+n-1),zpos(p+n-1),ox(1,m1),oy(1,m1),oz(1,m1),tx,ty,tz)
	    ! 'Zero' its position with respect to the axis centre...
	    tx=tx-ox(sp1,m1)
	    ty=ty-oy(sp1,m1)
	    tz=tz-oz(sp1,m1)
	    ! Transform the coordinate into the local coordinate system...
	    px=tx*axisx(sp1,m1,1) + ty*axisx(sp1,m1,2) + tz*axisx(sp1,m1,3)
	    py=tx*axisy(sp1,m1,1) + ty*axisy(sp1,m1,2) + tz*axisy(sp1,m1,3)
	    pz=tx*axisz(sp1,m1,1) + ty*axisz(sp1,m1,2) + tz*axisz(sp1,m1,3)
	    ! Accumulate the position....
	    avcat(n,1)=avcat(n,1)+px
	    avcat(n,2)=avcat(n,2)+py
	    avcat(n,3)=avcat(n,3)+pz
	    ! Intramolecular distribution.....
	    n1=NINT(px/0.15)
	    n2=NINT(py/0.15)
	    n3=NINT(pz/0.15)
	    ication(n1,n2,n3)=ication(n1,n2,n3)+1
	  end do
	  p=p+s_natoms(sp)
	end do
	! Next frame
	GOTO 101

120	write(0,*) "Finished."
	! * Normalise the data *
	! **********************
	do l=1,1	! Should be 1,nspecies but don't use anion-X data.
	  do k=1,nspecies
	write(0,*) "PDENSN",l,k,pdensn(l,k)
	pdensn(l,k) = pdensn(l,k) / nframes / s_nmols(k)
	write(0,*) "  (Reduced = ",pdensn(l,k),")"
	   do n1=-GRID,GRID
	    do n2=-GRID,GRID
	     do n3=-GRID,GRID
	      pdens(l,k,n1,n2,n3)=pdens(l,k,n1,n2,n3)/nframes/pdensn(l,k)/(delta**3)
	      ! pdens(l,k,n1,n2,n3)=pdens(l,k,n1,n2,n3)/nframes/200/(delta**3)
	     end do
	    end do
	   end do
	  end do
	end do

	do l=1,1
	  do k=1,MAXS
	    write(0,*) "Original pdensa,n",l,k,pdensa(l,k),pdensn(l,k)
	write(0,*) s_nmols,nframes
	    pdensa(l,k)=pdensa(l,k)/nframes/s_nmols(k)/s_nmols(l)
	  end do
	end do
	do n=1,s_natoms(1)
	  do m=1,3
	    avcat(n,m)=avcat(n,m)/nframes/s_nmols(1)
	  end do
	end do
	GOTO 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	GOTO 999

798	write(0,*) "Problem with OUTPUT file."
	GOTO 999
799	write(0,*) "HISTORY file ended prematurely!"
	GOTO 801
800	write(0,*) "End of unformatted HISTORY file found."
801	write(0,*) "Averages taken over ",nframes, " frames and ",pdensn(1,1)," atoms/species."
	write(0,*) "Numerical density averages (caught within GRID boundaries):"
	do n=1,1
	  do m=1,MAXS
	    write(0,*) "Species ",n," - Species ",m,": ",real(pdensa(n,m))/real((GRID*2+1)**3)
	  end do
	end do
	write(0,*) "Grid = ",(GRID*2+1),"x",(GRID*2+1),"x",(GRID*2+1)
	write(0,*) ""

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  IF (infile(n:n).EQ.".") THEN
	    baselen=n
	    GOTO 802
	  endIF
	end do
802	IF (baselen.EQ.-1) THEN
	  basename="3ddistresults."
	  baselen=14
	else
	  basename=infile(1:baselen)
	endIF

	do l=1,1	! Should be 1,nspecies, but anion is species 2 and we don't need that data.
	  do k=1,MAXS
	    resfile=basename(1:baselen)//"3ddist"//CHAR(48+l)//CHAR(48+k)
	    OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	    do n=-GRID,GRID
	      do m=-GRID,GRID
		do o=-GRID,GRID
		  ! write(9,*) n*delta,m*delta,o*delta,pdens(l,k,n,m,o)
		  write(9,*) pdens(l,k,n,m,o)
		end do
	      end do
	    end do
	    CLOSE(9)
	  end do
	end do
	resfile=basename(1:baselen)//"avcat"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	! Write out a pdb of the average cation in a format that the CPMS module can understand
	do n=1,s_natoms(1)
900	FORMAT ('ATOM  ',I5,A4,'0','UNK','I',I4,'R',3x,3F8.3)
	  write(9,"('Atom',4x,I3,2x,A1,15x,3(F7.4,2x))") n,atmnam(s_start(1)-1+n)(1:1),avcat(n,1),avcat(n,2),avcat(n,3)
	!ATOM      1  C   ACh     1    5.4012   64.3027 61.5836
	end do
	CLOSE(9)
	resfile=basename(1:baselen)//"icat"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=-GRID,GRID
	  do m=-GRID,GRID
	    do o=-GRID,GRID
	      write(9,*) n*delta,m*delta,o*delta,ication(n,m,o)
	    end do
	  end do
	end do
	CLOSE(9)
	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end

	subroutine pbc(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	use dlprw; implicit none
	! Performs minimum image convention.
	! xyz1: point to consider; xyz2: reference point; xyz3: minimum image coords of xyz1 (result)
	real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3
	! Performs minimum image convention.....
	x3=x1 - cell(1)*NINT((x1-x2)/cell(1))
	y3=y1 - cell(5)*NINT((y1-y2)/cell(5))
	z3=z1 - cell(9)*NINT((z1-z2)/cell(9))
	end subroutine pbc

	subroutine calcgc
	! Calculates the geometric/mass centres of all molecules in all species
	! s: element of GCxyz array to fill; sp: species no; natms: no. of atoms in gcal()
	use dlprw; use pdens_data
	implicit none
	real*8 :: discardr,massnorm
	real*8 tx,ty,tz,x0,y0,z0
	integer sp,n,m,o,p,comtype
	do sp=1,nspecies
	  do n=1,s_nmols(sp)
	    gcx(sp,n)=0.0
	    gcy(sp,n)=0.0
	    gcz(sp,n)=0.0
	  end do
	  if (s_natoms(sp).GT.1) then   ! Calculate centre of mass
	    o=s_start(sp)
	    do n=1,s_nmols(sp)
	      ! We will use zero as the reference point....
	      x0=0.0
	      y0=0.0
	      z0=0.0
	      massnorm=0.0
	      do m=1,s_natoms(sp)
	        CALL PBC(x(o+m-1),y(o+m-1),z(o+m-1),gcx(sp,n),gcy(sp,n),gcz(sp,n),tx,ty,tz,cell)
	        x0=x0 + (tx-gcx(sp,n))*mass(o+m-1)
	        y0=y0 + (ty-gcy(sp,n))*mass(o+m-1)
	        z0=z0 + (tz-gcz(sp,n))*mass(o+m-1)
	        massnorm=massnorm+mass(o+m-1)
	      end do
	      ! Normalise...
	      gcx(sp,n)=gcx(sp,n) + x0/massnorm
	      gcy(sp,n)=gcy(sp,n) + y0/massnorm
	      gcz(sp,n)=gcz(sp,n) + z0/massnorm
	      o=o+s_natoms(sp)
	    end do
	  end if
	end do
	end subroutine calcgc

	subroutine genaxis
	use dlprw; use pdens_data
	implicit none
	! Subroutine to calculate individual molecular axes of species molecules
	real*8 :: vx,vy,vz,tx,ty,tz
	real*8 :: vmag
	integer :: sp,o,n
	do sp=1,nspecies
	  if (s_natoms(sp).GT.2) THEN   ! Calculate local axes if we can!
	    o=s_start(sp)
	      do n=1,s_nmols(sp)
	      ! 1) Determine the X axis of the molecule
	      !    -- Get the vector between atoms 1 and 2...
	      call pbc(x(o+aa(sp,2)-1),y(o+aa(sp,2)-1),z(o+aa(sp,2)-1), &
		& x(o+aa(sp,1)-1),y(o+aa(sp,1)-1),z(o+aa(sp,1)-1),tx,ty,tz)
	      call getvector(tx,ty,tz,x(o+aa(sp,1)-1),y(o+aa(sp,1)-1),z(o+aa(sp,1)-1),vx,vy,vz)
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
	      call pbc(x(o+aa(sp,3)-1),y(o+aa(sp,3)-1),z(o+aa(sp,3)-1),ox(sp,n),oy(sp,n), &
		& oz(sp,n),tx,ty,tz)
	      call getvector(ox(sp,n),oy(sp,n),oz(sp,n),tx,ty,tz,vx,vy,vz)
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
	      !write(88,*) "O: ",ox(sp,n),oy(sp,n),oz(sp,n)
	      !write(88,*) "X: ",axisx(sp,n,1),axisx(sp,n,2),axisx(sp,n,3)
	      !write(88,*) "Y: ",axisy(sp,n,1),axisy(sp,n,2),axisy(sp,n,3)
	      !write(88,*) "Z: ",axisz(sp,n,1),axisz(sp,n,2),axisz(sp,n,3)
	    end do
	  else
	    ! Otherwise, set the origin of the species to the coords of the first atom.
	    ! For, e.g., halide ions, this is a must!
	    o=s_start(sp)
	    do n=1,s_nmols(sp)
	      ox(sp,n)=x(o)
	      oy(sp,n)=y(o)
	      oz(sp,n)=z(o)
	      o=o+s_natoms(sp)
	    end do
	  end if
	end do
	end subroutine genaxis

	subroutine getvector(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	implicit none
	real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3
	! xyz1,xyz2: coordinates of two atoms; xyz3: vector between them (result)
	real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3
	x3=x1-x2
	y3=y1-y2
	z3=z1-z2
	end subroutine getvector

	subroutine alloc_vars
	use pdens_data; use dlprw
	implicit none
	integer :: i
	! Geometric COM arrays
	i = max(s_nmols)
	allocate(gcx(nspecies,i)); allocate(gcy(nspecies,i)); allocate(gcz(nspecies,i))
	! Axis coordinate arrays
	allocate(axisx(nspecies,i,3)); allocace(axisy(nspecies,i,3)); allocace(axisz(nspecies,i,3))
	allocate(axisox(nspecies,i)); allocate(axisoy(nspecies,i)); allocate(axisoz(nspecies,i))

	end subroutine alloc_vars
