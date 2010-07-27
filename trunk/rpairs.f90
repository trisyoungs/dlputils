!	** rpairs.f90 **
!	Program to calculate the 3-dimensional distribution of different species about
!	one species' centres, and determine unique pair geometries based on local environment

	program rpairs
	use dlprw; use utility
	implicit none
	real*8, allocatable :: pdens(:,:,:)
	real*8, allocatable :: avggeom(:,:)		! Average species coordinates
	real*8, allocatable :: rrot(:,:)		! Rotated coordinates of axis centres
	integer, allocatable :: pdensi(:,:,:)		! Flag voxel array
	integer :: nfound, ncaught			! Total and 'binned' mols
	integer, allocatable :: molflags(:)		! Per-frame map of central species mols to include in averaging
	integer :: spexp				! Expected count
	integer :: grid = 20				! Grid in each direction for 3d grid
	integer :: minlobesize = 8			! Minimum lobesize (in voxels)
	real*8 :: delta = 0.5				! Default grid spacing
	logical :: molmap = .FALSE.			! Whether we're using a file of mol flags
	integer :: nmapframes, baselen,nargs,nframesused
	character*80 :: hisfile,outfile,basename,resfile,temp,flagfile,altheaderfile
	logical :: altheader = .FALSE.
	integer :: success,n1,n2,n3, n, nframes, m, o, p
	real*8 :: px, py, pz, mindist, maxdist, dist, tx, ty, tz, maximum, totalint, cutoff
	integer :: centresp,outersp,sp,m1,m2,startf,endf,n3dcentres,nlobes,lobesize
	integer :: iargc, lobeselect

	write(0,*) "*** rpairs"

	nargs = iargc()
	if (nargs.lt.2) then
	  write(*,"(a)") "Usage: rpairs <HIStory file> <OUTput file> [-options]"
	  write(*,"(a)") "        [-centre sp]            Set species sp to be the central one"
	  write(*,"(a)") "        [-outer sp]             Set species sp to be the outer one"
	  write(*,"(a)") "        [-axis sp x1 x2 y1 y2]  Atoms to use for axis calculation in species sp"
	  write(*,"(a)") "        [-grid npoints]         Grid points in each direction for prob. densities (default = 20)"
	  write(*,"(a)") "        [-delta spacing]        Spacing between grid points (default = 0.5 Angstrom)"
	  write(*,"(a)") "        [-minlobesize nvoxels]  Minimum allowable lobe size"
	  write(*,"(a)") "        [-start frameno]        Trajectory frame to start calculations (default = 1)"
	  write(*,"(a)") "        [-end frameno]          Trajectory frame to end calculations on (default = last)"
	  write(*,"(a)") "        [-molmap <file>]        Formatted file (I5,1X,n(I2,1x) of frame,mols to use of centresp"
	  write(*,"(a)") "        [-mindist r]            Minimum separation allowed between molecules (default = 0.0)"
	  write(*,"(a)") "        [-maxdist r]            Maximum separation allowed between molecules (default = 10000.0)"
	  write(*,"(a)") "        [-header file]          Use specified file to get header"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
	end if
	n = 10
	write(0,"(A,A)") "History file : ",hisfile
	write(0,"(A,A)") " Output file : ",outfile
	if (outinfo(outfile,1).eq.-1) goto 798

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  IF (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 50
	  endIF
	end do
50	if (baselen.eq.-1) then
	  basename="3ddistresults."
	  baselen=14
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=15,file=basename(1:baselen)//"out",form="formatted",status="replace")

	! Set some variable defaults before we read in any command-line arguments
	call alloc_axis()
	mindist = 0.0
	maxdist = 10000.0
	centresp = 1
	outersp = 2
	startf = 1
	endf = 0
	n3dcentres = 0

	n = 2
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-axis") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp,2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp,3)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp,4)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp," calculated from: X=",aa(sp,1),"->", &
	        & aa(sp,2),", Y=0.5X->0.5(r(",aa(sp,3),")->r(",aa(sp,4),"))"
	      axisdefined(sp) = .true.
	      if (aa(sp,1).eq.aa(sp,2)) then
		write(0,*) " --- X-axis atom ids are identical - atom will be used as centre but NO axes will be defined."
	        axisdefined(sp) = .false.
		axisusecom(sp) = .false.
	      end if
	    case ("-grid")
	      n = n + 1; call getarg(n,temp); read(temp,"(I4)") grid
	      write(0,"(A,I4)") "Grid points in each XYZ = ",grid
	    case ("-delta")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") delta
	      write(0,"(A,F6.3)") "Grid spacing = ",delta
	    case ("-minlobesize")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") minlobesize
	      write(0,"(A,I5)") "Minimum lobe size = ",minlobesize
	    case ("-centre")
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") centresp
	      write(0,"(A,I4)") "Central species is = ",centresp
	    case ("-outer")
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") outersp
	      write(0,"(A,I4)") "Outer species is = ",outersp
	    case ("-start")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") startf
	      write(0,"(A,I5)") "Starting frame = ",startf
	    case ("-end")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") endf
	      write(0,"(A,I5)") "End frame = ",endf
	    case ("-molmap")
	      n = n + 1; call getarg(n,flagfile)
	      write(0,"(A,A)") "Using flags for central species molecules from : ",flagfile
	      open(unit=20,file=flagfile,form="formatted",status="old")
	      molmap = .TRUE.
	    case ("-mindist")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") mindist
	      write(0,"(A,f8.4)") "Minimum separation between molecules = ",mindist
	      write(15,"(A,f8.4)") "Minimum separation between molecules = ",mindist
	    case ("-maxdist")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") maxdist
	      write(0,"(A,f8.4)") "Maximum separation between molecules = ",maxdist
	      write(15,"(A,f8.4)") "Maximum separation between molecules = ",maxdist
	    case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
              altheader = .TRUE.
	    case default
	      write(0,*) "Unrecognised command line option : ",temp
	      stop
	  end select
	end do
	
	! Print out a summary of the control variables to the output file.
	write(15,"(A,A)") "Input file: ",hisfile
	write(15,"(A,I3)") "Molecular species in file : ",nspecies
	do sp=1,nspecies
	  ! Check that molecules have had their axes defined
	  if (axisdefined(sp)) then
	    write(15,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Axis for species ",sp," calculated from : X=",aa(sp,1),"->", &
	      & aa(sp,2),", Y=0.5X->0.5(",aa(sp,3),"->",aa(sp,4),")"
	  else if (sp.eq.centresp) then
	    stop "A proper set of axes must be defined for the central species."
	  end if
	end do
	write(15,"(A,I4)") "Grid (intermolecular) points in each XYZ = ",grid
	write(15,"(A,F6.3)") "Grid spacing = ",delta
	write(15,"(A,I4)") "Central species = ",centresp
	write(15,"(A,I4)") "Outer species = ",outersp
	write(15,"(A,I5,A,I5,A)") "Frame range = ",startf," to ",endf," (0=last)"

	! Probability density arrays
	allocate(pdens(-grid:grid,-grid:grid,-grid:grid))
	allocate(pdensi(-grid:grid,-grid:grid,-grid:grid))
	! Rotated coordinates
	allocate(rrot(maxmols,3))
	! Average species etc.
	allocate(avggeom(maxatoms,3))
	if (molmap) allocate(molflags(s_nmols(centresp)))

	! Clear the arrays before we start...
	write(0,*) "Clearing initial arrays..."
	pdens = 0.0
	pdensi = 0
	ncaught = 0
	nfound = 0
	avggeom = 0.0

	! Read the header of the history file...
        call openhis(hisfile,10)
        if (readheader().EQ.-1) then
          if (altheader) then
            write(0,*) "Restarted trajectory:"
            close(dlpun_his)
            call openhis(altheaderfile,10)
            if (readheader().EQ.-1) goto 799
            close(dlpun_his)
            call openhis(hisfile,10)
          else
	    rewind(dlpun_his)
            !goto 799
          end if
        end if

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
	nframes = 0
	nframesused = 0
	nmapframes = 0
	! If we're using a mapfile, read it here
101	if (molmap) then
	  read(20,"(I5,1x,40(I2,1x))",end=118,err=117) m1,(molflags(m2),m2=1,s_nmols(centresp))
	  ! if (m1.NE.nframes) stop "Mapfile missed a frame number!"
	end if
102	success=readframe()
	IF (success.eq.1) goto 120  ! End of file encountered....
	IF (success.lt.0) goto 119  ! File error....
	nframes=nframes+1

	if (MOD(nframes,100).eq.0) then
	  if (molmap) then
	    write(0,"(I6,2x,'(',I6,',',I10,')')") nframes,nmapframes,n3dcentres
	  else
	    write(0,"(i6,2x,'(',I10,')')") nframes,n3dcentres
	  end if
	end if
	if (nframes.LT.startf) goto 102

	if ((molmap).and.(m1.ne.nframes)) goto 102
	if (molmap) nmapframes = nmapframes + 1

	nframesused = nframesused + 1

	! Generate all molecular axes
	call genaxis

	!
	! Calculate 3D distribution of outer species around central species
	!
	do m1=1,s_nmols(centresp)      ! Loop over all molecules of species 1...
	  ! If we're using a mapfile, decide whether to include this molecule
	  if (molmap) then
	    if (molflags(m1).eq.0) cycle 
	  end if
	  n3dcentres = n3dcentres + 1	! Normalisation counter
	  do m2=1,s_nmols(outersp)

	    ! Calculate minimum image position of molecule m2 about 0,0,0
	    call pbc(axisox(outersp,m2),axisoy(outersp,m2),axisoz(outersp,m2),axisox(centresp,m1),axisoy(centresp,m1),axisoz(centresp,m1),tx,ty,tz)
	    px = tx - axisox(centresp,m1)
	    py = ty - axisoy(centresp,m1)
	    pz = tz - axisoz(centresp,m1)
	    tx = px*axisx(centresp,m1,1) + py*axisx(centresp,m1,2) + pz*axisx(centresp,m1,3)
	    ty = px*axisy(centresp,m1,1) + py*axisy(centresp,m1,2) + pz*axisy(centresp,m1,3)
	    tz = px*axisz(centresp,m1,1) + py*axisz(centresp,m1,2) + pz*axisz(centresp,m1,3)

	    ! Check minimum / maximum separation
	    dist = sqrt(tx*tx + ty*ty + tz*tz)
	    if (dist.lt.mindist) cycle
	    if (dist.gt.maxdist) cycle

	    ! Calculate integer position in grid
	    n1=NINT(tx/delta)
	    n2=NINT(ty/delta)
	    n3=NINT(tz/delta)

	    ! If any of the n's are over grid, then only increase the counter
	    if (MAX(ABS(n1),MAX(ABS(n2),ABS(n3))).GT.grid) then
	      nfound = nfound + 1      ! No position, just count it....
	    else
	      if ((centresp.eq.outersp).and.(m1.eq.m2)) then
		! Same molecule, so don't add at all
	      else
	        ! Check to make sure the same molecule isn't being consider with itself  XXXXX
	        pdens(n1,n2,n3) = pdens(n1,n2,n3) + 1
	        nfound = nfound + 1
	        ncaught = ncaught + 1
	      end if
	    end if
	    p=p+s_natoms(outersp)
	  end do
	end do

	! Accumulate average geometry for central molecule
	p=s_start(centresp)
	do m1=1,s_nmols(centresp)
	  ! If we're using a mapfile, decide whether to include this molecule
	  if (molmap) then
	    if (molflags(m1).eq.0) then
	      p = p+s_natoms(centresp)
	      cycle
	    end if
	  end if
	  do n=1,s_natoms(centresp)
	    ! PBC the atom
	    call pbc(xpos(p+n-1),ypos(p+n-1),zpos(p+n-1),axisox(centresp,m1), &
		& axisoy(centresp,m1),axisoz(centresp,m1),tx,ty,tz)
	    ! 'Zero' its position with respect to the axis centre...
	    tx=tx-axisox(centresp,m1)
	    ty=ty-axisoy(centresp,m1)
	    tz=tz-axisoz(centresp,m1)
	    ! Transform the coordinate into the local coordinate system...
	    px=tx*axisx(centresp,m1,1) + ty*axisx(centresp,m1,2) + tz*axisx(centresp,m1,3)
	    py=tx*axisy(centresp,m1,1) + ty*axisy(centresp,m1,2) + tz*axisy(centresp,m1,3)
	    pz=tx*axisz(centresp,m1,1) + ty*axisz(centresp,m1,2) + tz*axisz(centresp,m1,3)
	    ! Accumulate the position....
	    avggeom(n,1)=avggeom(n,1)+px
	    avggeom(n,2)=avggeom(n,2)+py
	    avggeom(n,3)=avggeom(n,3)+pz
	  end do
	  p=p+s_natoms(centresp)
	end do

	! Next frame (or finish)
116	if (nframes.EQ.endf) goto 120
	goto 101

117	write(0,*) "Error reading mapfile."
	write(0,*) "Selected ",nmapframes," from trajectory before error."
	goto 120
118	write(0,*) "Reached end of mapfile."
	write(0,*) "Selected ",nmapframes," from trajectory."
	goto 120
119	write(0,*) "HISTORY file ended prematurely!"
	write(0,"(A,I5,A)") "Managed ",nframesused," frames before error."
	write(15,"(A)") "HISTORY file ended prematurely!"
	write(15,"(A,I5,A)") "Managed ",nframesused," frames before error."
120	write(0,*) "Finished."
	write(15,"(A)") "Finished."
	write(15,"(A,I5,A)") "Averages will be taken over ",nframesused," frames."

	! ### Normalise the data

	write(0,"(a,i1,a,i7,a,i7,a)") "Central species consisted ",n3dcentres," molecules"

	! Normalisation over frames

	! Calculate expected species numbers
	spexp = (s_nmols(centresp) - 1) * s_nmols(centresp) * nframesused

	! Species density about central species
	do n1=-grid,grid
	  do n2=-grid,grid
	    do n3=-grid,grid
	      pdens(n1,n2,n3)=pdens(n1,n2,n3)/n3dcentres/(delta**3)
	    end do
	  end do
	end do
 
	! Average molecule
	avggeom = avggeom / n3dcentres

	goto 801

	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "Problem with HISTORY file - failed to read header."
	goto 999
	write(0,*) "End of unformatted HISTORY file found."
	write(15,"(A)") "End of unformatted HISTORY file found."
801	write(0,"(A,I1)") " Central species = ",centresp
	write(0,"(A,I1,A,A)") " Species ",centresp,", ",s_name(centresp)
	write(0,"(A12,I9,A,I9,A)") "Expected : ",spexp," over ",nframesused," frames."
	write(0,"(A12,I9)") "Found : ",nfound
	write(0,"(A12,I9,A)") "Caught : ",ncaught," (in grid)"
	write(0,"(A12,I9)") "Selected : ",n3dcentres
	write(15,"(A,I1,A,A)") " Species ",centresp,", ",s_name(centresp)
	write(15,"(A12,I9,A,I9,A)") "Expected : ",spexp," over ",nframesused," frames."
	write(15,"(A12,I9)") "Found : ",nfound
	write(15,"(A12,I9,A)") "Caught : ",ncaught," (in grid)"
	write(15,"(A12,I9)") "Selected : ",n3dcentres
	write(0,"(3(A,I3),A,F4.1,A)") "Grid = ",(grid*2+1),"x",(grid*2+1),"x", &
	  & (grid*2+1)," ( x ",delta," Angstrom )"
	write(15,"(3(A,I3),A,F4.1,A)") "Grid = ",(grid*2+1),"x",(grid*2+1),"x", &
	  & (grid*2+1)," ( x ",delta," Angstrom )"
	write(0,"(A)") "-- Writing output files..."

	! -- Probability density of this molecule about the central species
	resfile=basename(1:baselen)//CHAR(48+centresp)//CHAR(48+sp)//".pdens"
	open(unit=9,file=resfile,form="formatted")
	write(9,*) 2*grid+1, 2*grid+1, 2*grid+1
	write(9,"(9f8.4)") delta,0.0,0.0,0.0,delta,0.0,0.0,0.0,delta
	write(9,"(3f10.4)") -grid*delta,-grid*delta,-grid*delta
	write(9,*) "zyx"
	do n=-grid,grid
	  do m=-grid,grid
	    do o=-grid,grid
	      write(9,"(f12.8)") pdens(n,m,o)
	    end do
	  end do
	end do
	close(9)

	! -- Average molecule
	resfile=basename(1:baselen)//"avg"//char(48+centresp)//".xyz"
	open(unit=9,file=resfile,form="formatted")
	! Write out pdb files of the average molecules
900	FORMAT (a2,5x,3(F9.5,2x))
	write(9,*) s_natoms(centresp)
	write(9,*) "Average: ", s_name(centresp)
	do n=1,s_natoms(centresp)
	  write(9,900) atmname(s_start(centresp)-1+n)(1:2),(avggeom(n,m),m=1,3)
	end do
	close(9)

	!
	! Analyse calculated density at different cutoffs to establish lobes
	!

	! Determine maximum value in grid and total integral
	maximum = -10000.0
	totalint = 0.0
	do n=-grid,grid
	  do m=-grid,grid
	    do o=-grid,grid
	      totalint = totalint + pdens(n,m,o)
	      if (pdens(n,m,o).gt.maximum) maximum = pdens(n,m,o)
	    end do
	  end do
	end do
	totalint = totalint * (delta**3)
	write(0,"(a,e12.5)") "Maximum value found in data is : ", maximum
	write(0,"(a,f10.6)") "Cubic volume of grid voxel : ", delta**3
	write(0,"(a,f10.6)") "Total integral of grid : ", totalint
	write(15,"(a,e12.5)") "Maximum value found in data is : ", maximum
	write(15,"(a,f10.6)") "Cubic volume of grid voxel : ", delta**3
	write(15,"(a,f10.6)") "Total integral of grid : ", totalint

	cutoff = maximum/2.0

	!do
	!  if (cutoff.gt.maximum) exit
	  ! Create subdensity of main grid, containing only values >= cutoff
	  pdensi = 0
	  do n=-grid,grid
	    do m=-grid,grid
	      do o=-grid,grid
	        if (pdens(n,m,o).ge.cutoff) pdensi(n,m,o) = -1
	      end do
	    end do
	  end do
	  write(0,*) sum(pdensi)
	  ! Determine how many distinct lobes there are in the quantised subgrid
	  nlobes = 0
	  ! Find an unselected grid voxel (== -1)
	  do n=-grid,grid
	    do m=-grid,grid
	      do o=-grid,grid
	        if (pdensi(n,m,o).ne.-1) cycle
		nlobes = nlobes + 1
		lobesize = lobeselect(pdensi,grid,n,m,o,nlobes)
	write(0,*) "LOBE",nlobes, lobesize
	      end do
	    end do
	  end do
	  
	  
	  cutoff = cutoff + maximum/20.0
	!end do

	! TEST - Write integer lobe map
	resfile=basename(1:baselen)//CHAR(48+centresp)//CHAR(48+sp)//"lobes.pdens"
	open(unit=9,file=resfile,form="formatted")
	write(9,*) 2*grid+1, 2*grid+1, 2*grid+1
	write(9,"(9f8.4)") delta,0.0,0.0,0.0,delta,0.0,0.0,0.0,delta
	write(9,"(3f10.4)") -grid*delta,-grid*delta,-grid*delta
	write(9,*) "zyx"
	do n=-grid,grid
	  do m=-grid,grid
	    do o=-grid,grid
	      write(9,"(f5.1)") pdensi(n,m,o)*1.0
	    end do
	  end do
	end do
	close(9)


	write(0,*) "Finished!"
	write(15,"(A)") "Finished!"
999	close(10)
	close(13)

	end program rpairs

	recursive integer function lobeselect(pdensi,grid,n,m,o,lobeid) result(res)
	implicit none
	integer :: grid, lobeid, n, m, o, res
	integer, intent(inout) :: pdensi(-grid:grid,-grid:grid,-grid:grid)
	integer :: nselected
	nselected = 0
	res = 0
	! Select target gridpoint
	if (pdensi(n,m,o).gt.0) then
	  write(0,*) "Oddness - gridpoint is already selected"
	  return
	end if
	pdensi(n,m,o) = lobeid
	nselected = nselected + 1
	! Check each neighbouring point in the grid - if it is -1, call the function again
	write(0,"(a,3i4,' to lobe ',i4)") "Added point ", n,m,o, lobeid
	if (n.gt.-grid) then
	  if (pdensi(n-1,m,o).eq.-1) nselected = nselected + lobeselect(pdensi,grid,n-1,m,o,lobeid)
	end if
	if (n.lt.grid) then
	  if (pdensi(n+1,m,o).eq.-1) nselected = nselected + lobeselect(pdensi,grid,n+1,m,o,lobeid)
	end if
	if (m.gt.-grid) then
	  if (pdensi(n,m-1,o).eq.-1) nselected = nselected + lobeselect(pdensi,grid,n,m-1,o,lobeid)
	end if
	if (m.lt.grid) then
	  if (pdensi(n,m+1,o).eq.-1) nselected = nselected + lobeselect(pdensi,grid,n,m+1,o,lobeid)
	end if
	if (o.gt.-grid) then
	  if (pdensi(n,m,o-1).eq.-1) nselected = nselected + lobeselect(pdensi,grid,n,m,o-1,lobeid)
	end if
	if (o.lt.grid) then
	  if (pdensi(n,m,o+1).eq.-1) nselected = nselected + lobeselect(pdensi,grid,n,m,o+1,lobeid)
	end if
	!write(0,*) nselected
	res = nselected
	end function lobeselect
