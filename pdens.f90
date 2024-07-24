!	** pdens.f90 **
!	Program to calculate the 3-dimensional distribution of different species about
!	one species' centres.
!	Updated version - rotates all molecules about target molecule first, then performs minimum
!	image before binning

	program pdensprog
	use dlprw; use utility; use parse; use PDensRW; use IList
	implicit none
	real*8, parameter :: radcon = 57.29577951d0
	! Target central and surrounding species
	integer :: centresp
	type(IntegerList) :: othersp
	! Intermolecular probability density
	type(PDens) :: pdensinter(MAXSP)		! Probability densities
	integer :: grid = 20				! Grid in each direction for intermolecular pdens
	real*8 :: delta = 0.5				! Default grid spacing for intermolecular pdens
	logical :: nointer = .FALSE.			! Whether intermolecular pdens calculation is enabled
	! Intramolecular probability density
	real*8, allocatable :: pdensintra(:,:,:)	! Intramolecular species distributions
	integer :: pgrid = 50				! Grid in each direction for intramolecular pdens
	real*8 :: pdelta = 0.15				! Default grid spacing for intramolecular pdens
	logical :: nointra = .TRUE.			! Whether intermolecular pdens calculation is enabled
	! Average geometry
	real*8, allocatable :: avggeom(:,:)		! Average species coordinates
	! Counters
	integer*8 :: npdenscentres, nintracentres
	integer*8 :: nfound(MAXSP), ncaught(MAXSP)	! Total and 'binned' mols
	integer*8 :: spexp(MAXSP)			! Expected species numbers
	! Map file for selecting central species
	logical :: molmap = .FALSE.			! Whether we're using a file of mol flags
	integer, allocatable :: molflags(:)		! Per-frame map of sp1 mols to include in averaging
	! Atomic sites for central species position
	type(IntegerList) :: originAtoms		! Atoms for central species origin (if any)
	! Atomic sites for surrounding species position
	type(IntegerList) :: atomSites(MAXSP)		! Atoms for surrounding species centres (if any)
	logical :: atomSitesAreCog(MAXSP) = .false.
	! Orientation checking
	logical :: orientcheck(MAXSP,3)			! Whether orientation check is enabled
	real*8 :: orientangle(20,3), orientdelta(20,3)	! Required orientations of molecules (if specified)
	! Selection, based on existing pdens
	logical :: selectsp(MAXSP)			! Whether selection based on old pdens is enabled
	real*8 :: selectmin(MAXSP)			! Minimum density value for a selection to succeed
	type(IntegerList) :: selectSites(MAXSP)		! Atom sites used to select positions in existing grid
	type(PDens) :: selectpdens(MAXSP)		! Reference pdens for 'select' mode
	character*80 :: selectfile(MAXSP)		! Filenames for old pdens
	! Frame control
	integer :: startf, endf

	! General variables
	integer :: nmapframes, baselen,nargs,nframesused
	character*80 :: hisfile,outfile,basename,resfile,temp,flagfile,altheaderfile
	logical :: altheader = .FALSE.
	integer :: success
	integer :: n, nframes, o, i, sp2index
	real*8 :: px, py, pz, mindistsq, maxdistsq, distsq, tx, ty, tz, v(3), origin(3)
	real*8 :: angdelta, ax
	integer :: sp1,m1,sp2,m2,p1,p2,nx,ny,nz,refnx,refny,refnz,maptot, mapint
	logical :: failed, getbin
	integer :: iargc

	write(0,*) "*** pdens"

	centresp = 0
	nargs = iargc()
	if (nargs.lt.3) then
	  write(*,"(a)") "Usage: pdens <HISTORYfile> <OUTPUTfile> <centresp> <othersp> [-options]"
	  write(*,"(a)") "        [-axis sp x1 x2 y1 y2]     Atoms to use for axis calculation in species sp"
	  write(*,"(a)") "        [-atoms sp <list>]         Add list of atoms in species sp to use as other points"
	  write(*,"(a)") "        [-cartesian sp]            Use Cartesian axes for species sp"
	  write(*,"(a)") "        [-cog sp <list>]           Use centre-of-geometry calculated from atom list for other point"
	  write(*,"(a)") "        [-delta spacing]           Spacing between grid points (default = 0.5 Angstrom)"
	  write(*,"(a)") "        [-end frameno]             Trajectory frame to end calculations on (default = last)"
	  write(*,"(a)") "        [-grid npoints]            Grid points in each direction for prob. densities (default = 20)"
	  write(*,"(a)") "        [-header file]             Use specified file to get header"
	  write(*,"(a)") "        [-maxdist r]               Maximum separation allowed between molecules (default = 10000.0)"
	  write(*,"(a)") "        [-mindist r]               Minimum separation allowed between molecules (default = 0.0)"
	  write(*,"(a)") "        [-molmap <file>]           Formatted file (I5,1X,n(I2,1x) of frame,mols to use of sp1"
	  write(*,"(a)") "        [-nointer]                 Prohibit calculation of the intermolecular probability densities"
	  write(*,"(a)") "        [-intra]                   Enable calculation of the intramolecular probability density and average molecule"
	  write(*,"(a)") "        [-orient sp axis angle delta]  Specify a restriction in the angular orientation of the second &
		& molecule (use negative delta for symmetry about 90deg)"
	  write(*,"(a)") "        [-origin <list>]           Use centre-of-geometry calculated from atom list for central point"
	  write(*,"(a)") "        [-otheraxis sp x1 x2 y1 y2] Alternative axis definition, for surrounding &
		& molecule position and orientation"
	  write(*,"(a)") "        [-pdelta spacing]          Spacing between pgrid points (default = 0.15 Angstrom)"
	  write(*,"(a)") "        [-pgrid npoints]           Grid points in each direction for species densities (default = 50)"
	  write(*,"(a)") "        [-select sp pdens atoms cut] Select only other centres for which the centre of geometry of the atoms supplied&
		& exists in a position in the supplied pdens with value greater than cutoff specified"
	  write(*,"(a)") "        [-start frameno]           Trajectory frame to start calculations (default = 1)"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
	  call getarg(3,temp); read (temp,"(i4)") centresp
	  call getarg(4,temp); if (.not.parseIntegerList(temp, othersp)) stop "Failed to parse othersp list."
	end if

	write(0,"(A,A)") "History file : ",hisfile
	write(0,"(A,A)") " Output file : ",outfile
	if (outinfo(outfile,1).eq.-1) goto 798
	write(0,"(A,i3)") "Central species is = ",centresp
	write(0,"(A,20i3)") "Other species are = ", othersp%items(1:othersp%n)

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

	! Set some variable defaults before we read in any command-line arguments
	call alloc_axis()
	mindistsq = 0.0
	maxdistsq = 1.0e9
	startf = 1
	endf = 0
	npdenscentres = 0
	nintracentres = 0
	orientcheck = .FALSE.
	orientangle = 0.0
	orientdelta = 0.0
	selectsp = .false.

	n = 4
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-atoms")
	      n = n + 1; sp1 = getargi(n)
	      if (atomSites(sp1)%n.gt.0) stop "Definition of sites (-atoms or -cog) specified twice for species"
	      n = n + 1; call getarg(n,temp); if (.not.parseIntegerList(temp, atomSites(sp1))) stop "Failed to parse atom list."
	      write(0,"(i2,a,i2,a,20i4)") atomSites(sp1)%n, " atoms added as other sites for species ", sp1, ": ", atomSites(sp1)%items(1:atomSites(sp1)%n)
	    case ("-axis") 
	      n = n + 1; sp1 = getargi(n)
	      n = n + 1; axesAatoms(sp1,1) = getargi(n)
	      n = n + 1; axesAatoms(sp1,2) = getargi(n)
	      n = n + 1; axesAatoms(sp1,3) = getargi(n)
	      n = n + 1; axesAatoms(sp1,4) = getargi(n)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp1," calculated from: X=",axesAatoms(sp1,1),"->", &
	        & axesAatoms(sp1,2),", Y=0.5X->0.5(r(",axesAatoms(sp1,3),")->r(",axesAatoms(sp1,4),"))"
	      do i=1,3
		if ((axesAatoms(sp1,i).lt.1).or.(axesAatoms(sp1,i).gt.s_natoms(sp1))) stop "Atom id out of range for axes on this species!"
	      end do
	      axesAdefined(sp1) = .true.
	    case ("-cartesian")
	      n = n + 1; sp1 = getargi(n)
	      axesAatoms(sp1,:) = -1
	      axesAdefined(sp1) = .true.
	      write(0,"(A,i3)") "Cartesian axes will be used for species ", sp1
	    case ("-cog")
	      n = n + 1; sp1 = getargi(n)
	      if (atomSites(sp1)%n.gt.0) stop "Definition of sites (-atoms or -cog) specified twice for species"
	      atomSitesAreCog(sp1) = .true.
	      n = n + 1; call getarg(n,temp); if (.not.parseIntegerList(temp, atomSites(sp1))) stop "Failed to parse atom list."
	      write(0,"(i2,a,i2,a,20i4)") atomSites(sp1)%n, " atoms added as COG site for species ", sp1, ": ", atomSites(sp1)%items(1:atomSites(sp1)%n)
	    case ("-delta")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") delta
	      write(0,"(A,F6.3)") "Grid spacing = ",delta
	    case ("-end")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") endf
	      write(0,"(A,I5)") "End frame = ",endf
	    case ("-grid")
	      n = n + 1; grid = getargi(n)
	      write(0,"(A,I4)") "Grid points in each XYZ = ",grid
	    case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
              altheader = .TRUE.
	    case ("-intra")
	      nointra = .false.
	      write(0,"(A,I4)") "Intramolecular probability distributions will be computed."
	    case ("-mindist")
	      n = n + 1; mindistsq = getargr(n)
	      write(0,"(A,f8.4)") "Minimum separation between molecules = ",mindistsq
	      mindistsq = mindistsq * mindistsq
	    case ("-maxdist")
	      n = n + 1; maxdistsq = getargr(n)
	      write(0,"(A,f8.4)") "Maximum separation between molecules = ",maxdistsq
	      maxdistsq = maxdistsq * maxdistsq
	    case ("-molmap")
	      n = n + 1; call getarg(n,flagfile)
	      write(0,"(A,A)") "Using flags for species 1 molecules from : ",flagfile
	      open(unit=20,file=flagfile,form="formatted",status="old")
	      molmap = .TRUE.
	    case ("-nointer")
	      nointer = .TRUE.
	      write(0,"(A,I4)") "Intermolecular probability distributions will not be computed."
	    case ("-orient")
	      n = n + 1; sp1 = getargi(n)
	      n = n + 1; i = getargi(n)
	      n = n + 1; orientangle(sp1,i) = getargr(n)
	      n = n + 1; orientdelta(sp1,i) = getargr(n)
	      write(0,"(A,i2,a,i1,a,f10.4,a,f10.4,a)") "Only molecules of sp ",sp1," with axis ", i, " angle delta of ", &
		& orientdelta(sp1,i), " degrees about ", orientangle(sp1,i), " will be binned"
	      orientcheck(sp1,i) = .TRUE.
	    case ("-origin")
	      n = n + 1; call getarg(n,temp); if (.not.parseIntegerList(temp, originAtoms)) stop "Failed to parse atom list for origin."
	      write(0,"(i2,a,20i4)") originAtoms%n, " atoms added as origin site for central species :", originAtoms%items(1:originAtoms%n)
	    case ("-otheraxis") 
	      n = n + 1; sp1 = getargi(n)
	      if (atomSites(sp1)%n.ne.0) stop "Definition of sites and otheraxis cannot be used together."
	      n = n + 1; axesBatoms(sp1,1) = getargi(n)
	      n = n + 1; axesBatoms(sp1,2) = getargi(n)
	      n = n + 1; axesBatoms(sp1,3) = getargi(n)
	      n = n + 1; axesBatoms(sp1,4) = getargi(n)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Alternative axes for species ",sp1," calculated from: X=",axesBatoms(sp1,1),"->", &
	        & axesBatoms(sp1,2),", Y=0.5X->0.5(r(",axesBatoms(sp1,3),")->r(",axesBatoms(sp1,4),"))"
	      axesBdefined(sp1) = .true.
	    case ("-pdelta")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") pdelta
	      write(0,"(A,F6.3)") "PGrid spacing = ",pdelta
	    case ("-pgrid")
	      n = n + 1; pgrid = getargi(n)
	      write(0,"(A,I4)") "PGrid points in each XYZ = ",pgrid
	    case ("-select")
	      n = n + 1; sp1 = getargi(n)
	      selectsp(sp1) = .true.
	      n = n + 1; call getarg(n,selectfile(sp1))
	      n = n + 1; call getarg(n,temp); if (.not.parseIntegerList(temp, selectSites(sp1))) stop "Failed to parse atom list for selection site."
	      n = n + 1; selectmin(sp1) = getargr(n)
	      write(0,"(a,i2,a,f6.2,a,20i4)") "Species ", sp1, "set for selection, cutoff = ", selectmin(sp1), "and based on COG of atoms ", selectSites(sp1)%items(1:selectSites(sp1)%n)
	    case ("-start")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") startf
	      write(0,"(A,I5)") "Starting frame = ",startf
	    case default
	      write(0,*) "Unrecognised command line option : ",temp
	      stop
	  end select
	end do
	
	! Print out a summary of the control variables to the output file.
	write(6,"(A,A)") "Input file: ",hisfile
	write(6,"(A,I3)") "Molecular species in file : ",nspecies
	do sp1=1,nspecies
	  ! Check that the necessary molecules have had their axes defined
	  if (axesAdefined(sp1)) then
	    write(6,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp1," calculated from: X=",axesAatoms(sp1,1),"->", &
	      & axesAatoms(sp1,2),", Y=0.5X->0.5(r(",axesAatoms(sp1,3),")->r(",axesAatoms(sp1,4),"))"
	  else if (sp1.eq.centresp) then
	    stop "A proper set of axes must be defined for the central species."
	  else if (orientcheck(sp1,1).or.orientcheck(sp1,2).or.orientcheck(sp1,3)) then
	    stop "Axes must be defined on any species whose orientation is to be considered."
	  end if
	end do
	write(6,"(A,I4)") "Grid (intermolecular) points in each XYZ = ",grid
	write(6,"(A,I4)") "PGrid (intramolecular) points in each XYZ = ",pgrid
	write(6,"(A,F6.3)") "Grid spacing = ",delta
	write(6,"(A,F6.3)") "PGrid spacing = ",pdelta
	write(6,"(A,I4)") "Central species = ",centresp
	write(6,"(A,I5,A,I5,A)") "Frame range = ",startf," to ",endf," (0=last)"
	if (nointra.and.nointer) stop "Nothing to do (both -nointra and -nointer given)."

	! Probability density arrays
	allocate(pdensintra(-pgrid:pgrid,-pgrid:pgrid,-pgrid:pgrid))
	! Average species etc.
	allocate(avggeom(maxatoms,3))
	if (molmap) allocate(molflags(s_nmols(centresp)))

	! Setup arrays etc
	write(0,*) "Setting up initial arrays..."
	do n=1,othersp%n
	  sp1 = othersp%items(n)
	  pdensinter(sp1)%axes = (/ delta,0.0d0,0.0d0,0.0d0,delta,0.0d0,0.0d0,0.0d0,delta /)
	  pdensinter(sp1)%origin = (/ -grid*delta,-grid*delta,-grid*delta /)
	  if (.not.allocPDens(pdensinter(sp1),-grid,-grid,-grid,grid,grid,grid)) stop
	end do
	ncaught = 0
	nfound = 0
	pdensintra = 0.0
	spexp = 0
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
            goto 799
          end if
        end if

	! Load in and process any selection pdens files...
	do n=1,nspecies
	  if (.not.selectsp(n)) cycle
	  if (.not.loadPDens(selectfile(n), selectpdens(n))) stop "Failed to load reference pdens."
	  ! Make some basic checks on grid extent etc...
	  do i=1,3
	    if ((selectpdens(n)%gridMin(i).ne.pdensinter(centresp)%gridMin(1)).or.(selectpdens(n)%gridMax(i).ne.pdensinter(centresp)%gridMax(i))) then
	      write(0,*) "Error: Supplied reference pdens has different grid to that requested for this calculation."
	      write(0,"(a,6i3)") "Reference: ", selectpdens(n)%gridMin, selectpdens(n)%gridMax
	      write(0,"(a,6i3)") "Current: ", pdensinter(centresp)%gridMin,pdensinter(centresp)%gridMax
	      stop
	    end if
	  end do
	end do

	! XXXX
	! XXXX Main routine....
	! XXXX
	nframes = 0
	nframesused = 0
	nmapframes = 0
        sp1=centresp
	!if (molmap) molflags(:) = 0
	! If we're using a mapfile, read it here
        ! TFH EDITS HERE - YOU HAVE BEEN WARNED
        ! reads molmap file line until a "-1" is reached 
101	if (molmap) then
          molflags(:) = 0
	  read(20,"(I5)",advance="no",end=118,err=117) m1 !,(molflags(m2),m2=1,s_nmols(centresp))
	  if (m1.NE.(nframes+1)) stop "Mapfile missed a frame number!"
          do m2=1,s_nmols(centresp)
             read(20,"(1x,I5)",advance="no") mapint  !molflags(m2)
             !write(0,*) mapint
 	     if (mapint.GT.0) then
	        molflags(mapint)=1
             elseif (mapint.EQ.-1) then 
               exit
             else
               stop "Error in mapfile"
             endif               
          enddo
          read(20,*)
          !write(0,*) "for frame ,", m1, " molflags = ", molflags(:)!,m2=1,s_nmols(centresp))
	end if
102	success=readframe()
	IF (success.eq.1) goto 120  ! End of file encountered....
	IF (success.lt.0) goto 119  ! File error....
	nframes=nframes+1

	if (nframes.LT.startf) goto 102

	if ((molmap).and.(m1.ne.nframes)) goto 102
	if (molmap) nmapframes = nmapframes + 1

	nframesused = nframesused + 1

	! Generate all molecular axes
	call genaxes()

	!
	! Calculate 3D distributions
	!
	if (nointer) goto 105
	p1 = s_start(sp1) - s_natoms(sp1) - 1
	! Loop over all molecules of species 1...
	do m1=1,s_nmols(sp1)

	  p1 = p1 + s_natoms(sp1)

	  ! If we're using a mapfile, decide whether to include this molecule
	  if (molmap) then
	    if (molflags(m1).eq.0) cycle 
            !write(0,*) 'USing molecule ', m1, 'molflag = ', molflags(m1)
	  end if
	  npdenscentres = npdenscentres + 1	! Normalisation counter

	  ! Get origin on central molecule for distance measurements
	  if (originAtoms%n.gt.0) then
	    ! Get average coordinate, PBC it with the axis origin, and subtract the axis origin
	    call averagePosition(originAtoms%items,originAtoms%n,p1,v)
	    call pbc(v(1), v(2), v(3), axesAorigin(sp1,m1,1), axesAorigin(sp1,m1,2), axesAorigin(sp1,m1,3),tx,ty,tz)
	    tx = tx - axesAorigin(sp1,m1,1)
	    ty = ty - axesAorigin(sp1,m1,2)
	    tz = tz - axesAorigin(sp1,m1,3)

	    ! Rotate into reference frame
	    origin(1) = tx*axesA(sp1,m1,1) + ty*axesA(sp1,m1,2) + tz*axesA(sp1,m1,3)
	    origin(2) = tx*axesA(sp1,m1,4) + ty*axesA(sp1,m1,5) + tz*axesA(sp1,m1,6)
	    origin(3) = tx*axesA(sp1,m1,7) + ty*axesA(sp1,m1,8) + tz*axesA(sp1,m1,9)
	  else
	    origin = 0.0
	  end if

	  ! Loop over second (surrounding) species
	  do sp2index=1,othersp%n

	    sp2 = othersp%items(sp2index)

	    p2 = s_start(sp2) - s_natoms(sp2) - 1
	    do m2=1,s_nmols(sp2)

	      ! Increase atom offset
	      p2 = p2 + s_natoms(sp2)

	      ! Don't consider central molecule with itself
	      if ((sp1.eq.sp2).and.(m1.eq.m2)) cycle

	      ! Check orientation if it has been specified
	      failed = .false.
	      do n=1,3
	        if (orientcheck(sp2,n)) then
	          ! Take dot product of axis 'n' on sp2 with that of the central molecule, and calculate the angle delta
		  if (axesBdefined(sp2)) then
		    ax = axesB(sp2,m2,(n-1)*3+1)*axesA(sp1,m1,(n-1)*3+1) + axesB(sp2,m2,(n-1)*3+2)*axesA(sp1,m1,(n-1)*3+2) + axesB(sp2,m2,(n-1)*3+3)*axesA(sp1,m1,(n-1)*3+3)
		  else
		    ax = axesA(sp2,m2,(n-1)*3+1)*axesA(sp1,m1,(n-1)*3+1) + axesA(sp2,m2,(n-1)*3+2)*axesA(sp1,m1,(n-1)*3+2) + axesA(sp2,m2,(n-1)*3+3)*axesA(sp1,m1,(n-1)*3+3)
		  end if
	          angdelta = safeAngle(ax)
	          ! If delta is negative, map delta onto +/-90
	          if ((orientdelta(sp2,n).lt.0.0).and.(angdelta.gt.90.0)) angdelta = angdelta - 180.0
	          if (dabs(angdelta-orientangle(sp2,n)).gt.dabs(orientdelta(sp2,n))) failed = .true.
		end if
		if (failed) exit
	      end do
	      if (failed) cycle

	      ! Filter based on position/density of site in reference pdens?
	      if (selectsp(sp2)) then
		! Get target bin from atom list supplied
		call averagePosition(selectSites(sp2)%items,selectSites(sp2)%n,p2,v)
		if (getbin(sp1,m1,origin,0.0d0,maxdistsq*2.0,v(1),v(2),v(3),refnx,refny,refnz,grid,delta)) then
		  if (selectpdens(sp2)%grid(refnx,refny,refnz).lt.selectmin(sp2)) cycle
		end if
	      end if

	      ! Pass position of axesAorigin, or axesBorigin, or (cog of) individual atomSites, depending on what has been defined
	      ! If atomSites() have been defined, these override everything else
	      if ((atomSites(sp2)%n.gt.0).and.atomSitesAreCog(sp2)) then
		! Use list of supplied atoms, forming a COG with them
		call averagePosition(atomSites(sp2)%items, atomSites(sp2)%n, p2, v)
		if (getbin(sp1,m1,origin,mindistsq,maxdistsq,v(1),v(2),v(3),nx,ny,nz,grid,delta)) then
		  pdensinter(sp2)%grid(nx,ny,nz) = pdensinter(sp2)%grid(nx,ny,nz) + 1
		  ncaught(sp2) = ncaught(sp2) + 1
		end if
		nfound(sp2) = nfound(sp2) + 1
	      else if (atomSites(sp2)%n.gt.0) then
		! Use list of supplied atoms, with each contributing a point to the pdens
		do n=1,atomSites(sp2)%n
		  i = p2+atomSites(sp2)%items(n)
		  if (getbin(sp1,m1,origin,mindistsq,maxdistsq,xpos(i),ypos(i),zpos(i),nx,ny,nz,grid,delta)) then
		    pdensinter(sp2)%grid(nx,ny,nz) = pdensinter(sp2)%grid(nx,ny,nz) + 1
		    ncaught(sp2) = ncaught(sp2) + 1
		  end if
		  nfound(sp2) = nfound(sp2) + atomSites(sp2)%n
		end do
	      else if (axesBdefined(sp2)) then
		! Use origin of alternative axes
		if (getbin(sp1,m1,origin,mindistsq,maxdistsq,axesBorigin(sp2,m2,1),axesBorigin(sp2,m2,2),axesBorigin(sp2,m2,3),nx,ny,nz,grid,delta)) then
		  pdensinter(sp2)%grid(nx,ny,nz) = pdensinter(sp2)%grid(nx,ny,nz) + 1
		  ncaught(sp2) = ncaught(sp2) + 1
		end if
		nfound(sp2) = nfound(sp2) + 1
	      else
		! Use origin of axes
		if (getbin(sp1,m1,origin,mindistsq,maxdistsq,axesAorigin(sp2,m2,1),axesAorigin(sp2,m2,2),axesAorigin(sp2,m2,3),nx,ny,nz,grid,delta)) then
		  pdensinter(sp2)%grid(nx,ny,nz) = pdensinter(sp2)%grid(nx,ny,nz) + 1
		  ncaught(sp2) = ncaught(sp2) + 1
		end if
		nfound(sp2) = nfound(sp2) + 1
	      end if

	    end do
	  end do
	end do

	! Calculate average molecule and intramolecular distribution
105	p1 = s_start(sp1) - s_natoms(sp1) - 1
	do m1=1,s_nmols(sp1)

	  p1 = p1 + s_natoms(sp1)

	  ! If we're using a mapfile, decide whether to include this molecule
	  if (molmap.and.(molflags(m1).eq.0)) cycle

	  nintracentres = nintracentres + 1

	  ! Calculate origin offset (if there is one)
	  origin = 0.0
	  if (originAtoms%n.gt.0) then
	    ! Get average coordinate, PBC it with the axis origin, and subtract the axis origin
	    call averagePosition(originAtoms%items,originAtoms%n,p1,v)
	    call pbc(v(1), v(2), v(3), axesAorigin(sp1,m1,1), axesAorigin(sp1,m1,2), axesAorigin(sp1,m1,3),tx,ty,tz)
	    tx = tx - axesAorigin(sp1,m1,1)
	    ty = ty - axesAorigin(sp1,m1,2)
	    tz = tz - axesAorigin(sp1,m1,3)

	    ! Rotate into reference frame
	    origin(1) = tx*axesA(sp1,m1,1) + ty*axesA(sp1,m1,2) + tz*axesA(sp1,m1,3)
	    origin(2) = tx*axesA(sp1,m1,4) + ty*axesA(sp1,m1,5) + tz*axesA(sp1,m1,6)
	    origin(3) = tx*axesA(sp1,m1,7) + ty*axesA(sp1,m1,8) + tz*axesA(sp1,m1,9)
            
           ! if (m1.EQ.1) then
		!write(*,*) origin
                !write(*,*) axesAorigin(sp1,m1,:)
	   ! end if

  
	  
	  end if

	  ! Loop over atoms
	  do n=1,s_natoms(sp1)

	    ! PBC the atom
	    call pbc(xpos(p1+n),ypos(p1+n),zpos(p1+n),axesAorigin(sp1,m1,1), &
	  	& axesAorigin(sp1,m1,2),axesAorigin(sp1,m1,3),tx,ty,tz)

	    ! 'Zero' its position with respect to the axis centre...
	    tx=tx-axesAorigin(sp1,m1,1)+origin(1)
	    ty=ty-axesAorigin(sp1,m1,2)+origin(2)
	    tz=tz-axesAorigin(sp1,m1,3)+origin(3)

	    ! Transform the coordinate into the reference frame
	    px=tx*axesA(sp1,m1,1) + ty*axesA(sp1,m1,2) + tz*axesA(sp1,m1,3)
	    py=tx*axesA(sp1,m1,4) + ty*axesA(sp1,m1,5) + tz*axesA(sp1,m1,6)
	    pz=tx*axesA(sp1,m1,7) + ty*axesA(sp1,m1,8) + tz*axesA(sp1,m1,9)

	    !Add on origin offset
            ! TFH EDIT 08/10/2018 to "fix" that average molecule orgin not in the correct position if using -origin
            px=px-origin(1)
            py=py-origin(2)
            pz=pz-origin(3)

	    ! Accumulate the position....
	    avggeom(n,1) = avggeom(n,1)+px
	    avggeom(n,2) = avggeom(n,2)+py
	    avggeom(n,3) = avggeom(n,3)+pz

	    ! Intramolecular distribution.....
	    if (.not.nointra) then
	      nx=NINT(px/pdelta)
	      ny=NINT(py/pdelta)
	      nz=NINT(pz/pdelta)
	      if (MAX(ABS(nx),MAX(ABS(ny),ABS(nz))).lt.pgrid) then
	        pdensintra(nx,ny,nz) = pdensintra(nx,ny,nz) + 1
	      else
	        write(0,*) "Large distance : nx,ny,nz = ",nx,ny,nz
	      end if
	    end if
	  end do

	end do

110	if (MOD(nframes,100).eq.0) then
	  if (molmap) then
	    write(0,"(I6,2x,'(',I6,',',I10,',',I10')')") nframes,nmapframes,npdenscentres,nintracentres
	  else
	    write(0,"(i6,2x,'(',I10,',',I10,')')") nframes,npdenscentres,nintracentres
	  end if
	end if

	! Next frame (or finish)
115	if (nframes.EQ.endf) goto 120
	goto 101

117	write(6,*) "Error reading mapfile."
	write(6,*) "Selected ",nmapframes," from trajectory before error."
	goto 120
118	write(6,*) "Reached end of mapfile."
	write(6,*) "Selected ",nmapframes," from trajectory."
	goto 120
119	write(6,*) "HISTORY file ended prematurely!"
	write(6,"(A,I5,A)") "Managed ",nframesused," frames before error."
120	write(6,*) "Finished."
	write(6,"(A,I5,A)") "Averages will be taken over ",nframesused," frames."

	! ### Normalise the data

	write(6,"(a,i7,a)") "Intermolecular calculation considers ",npdenscentres," molecules"
	write(6,"(a,i7,a)") "Intramolecular calculation considers ",nintracentres," molecules"

	! Set normalisation over frames and central species molecules
	do n=1,othersp%n

	  sp2 = othersp%items(n)

	  ! Calculate expected species numbers
	  spexp(sp2) = s_nmols(centresp) * s_nmols(sp2) * nframesused
	  if (centresp.eq.sp2) spexp(sp2) = spexp(sp2) - s_nmols(sp2) * nframesused

	  if ((atomSites(sp2)%n.gt.0).and..not.atomSitesAreCog(sp2)) spexp(sp2) = spexp(sp2) * atomSites(sp2)%n

	  ! Species density about central species
	  pdensinter(sp2)%grid(:,:,:) = pdensinter(sp2)%grid(:,:,:) / npdenscentres / (delta**3)

	end do
 
	! Intramolecular distribution
        pdensintra(:,:,:) = pdensintra(:,:,:)/nintracentres/(pdelta**3)

	! Average molecule
	avggeom(:,:) = avggeom(:,:)/nintracentres

	goto 801

	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "Problem with HISTORY file - failed to read header."
	goto 999
	write(6,"(A)") "End of unformatted HISTORY file found."
801	write(6,"(A,I1)") " Central species = ",centresp
	do sp1=1,nspecies
	  write(6,"(A,I1,A,A)") " Species ",sp1,", ", s_name(sp1)
	  write(6,"(A,I12,A,I9,A)") "Expected : ", spexp(sp1)," over ",nframesused," frames."
	  write(6,"(A,I12)") " Points tested : ", nfound(sp1)
	  write(6,"(A,f10.6)") " Points number density : ", nfound(sp1) / volume(cell)
	  write(6,"(A,I12,A)") "Caught : ", ncaught(sp1), " (in grid)"
	  if (sp1.eq.centresp) write(6,"(A,I9)") "Selected (inter) : ",npdenscentres
	  if (sp1.eq.centresp) write(6,"(A,I9)") "Selected (intra) : ",nintracentres
	end do
	write(6,"(3(A,I3),A,F4.1,A)") "Grid = ",(grid*2+1),"x",(grid*2+1),"x", &
	  & (grid*2+1)," ( x ",delta," Angstrom )"
	write(6,"(A)") "-- Writing output files..."

	do i=1,othersp%n

	  sp2 = othersp%items(i)

	  if (.not.nointer) then
	    ! -- Probability density of this molecule about the central species
	    resfile=basename(1:baselen)//CHAR(48+centresp)//CHAR(48+sp2)//".pdens"
	    if (.not.savePDens(resfile, pdensinter(sp2))) write(0,*) "ERROR: Failed to save data"
	  end if

	end do

	! -- Average molecule
	resfile=basename(1:baselen)//"avg"//char(48+centresp)//".xyz"
	open(unit=9,file=resfile,form="formatted")
	! Write out xyz files of the average molecules
900	FORMAT (a2,5x,3(F9.5,2x))
	write(9,*) s_natoms(centresp)
	write(9,*) "Average: ", s_name(centresp)
	do n=1,s_natoms(centresp)
	  write(9,900) atmname(s_start(centresp)-1+n)(1:2),(avggeom(n,i),i=1,3)
	end do
	close(9)

	if (.not.nointra) then
	  ! -- Intramolecular probability distribution
	  resfile=basename(1:baselen)//"intra"//char(48+centresp)//".pdens"
	  open(unit=9,file=resfile,form="formatted")
	  write(9,*) 2*pgrid+1, 2*pgrid+1, 2*pgrid+1
	  write(9,"(9f8.4)") pdelta,0.0,0.0,0.0,pdelta,0.0,0.0,0.0,pdelta
	  write(9,"(3f10.4)") -pgrid*pdelta,-pgrid*pdelta,-pgrid*pdelta
	  write(9,*) "zyx"
	  do nx=-pgrid,pgrid
	    do ny=-pgrid,pgrid
	      do nz=-pgrid,pgrid
	        write(9,"(F12.8)") pdensintra(nx,ny,nz)
	      end do
	    end do
	  end do
	  close(9)
	end if

	write(6,"(A)") "Finished!"
999	close(10)

	end program pdensprog

	logical function getbin(sp1,m1,origin,mindistsq,maxdistsq,x,y,z,nx,ny,nz,grid,griddelta)
	use utility
	implicit none
	integer, intent(in) :: sp1, m1, grid
	integer, intent(out) :: nx, ny, nz
	real*8, intent(in) :: x, y, z, mindistsq, maxdistsq, griddelta
	real*8 :: px, py, pz, tx, ty, tz, offx, offy, offz, origin(3), distsq

	! Get mim of supplied coordinate with the origin of the central species
	call pbc(x,y,z,axesAorigin(sp1,m1,1),axesAorigin(sp1,m1,2),axesAorigin(sp1,m1,3),tx,ty,tz)

	! Get vector between minimum image coordinates and central molecule axis origin
	px = tx - axesAorigin(sp1,m1,1)
	py = ty - axesAorigin(sp1,m1,2)
	pz = tz - axesAorigin(sp1,m1,3)
	tx = px*axesA(sp1,m1,1) + py*axesA(sp1,m1,2) + pz*axesA(sp1,m1,3)
	ty = px*axesA(sp1,m1,4) + py*axesA(sp1,m1,5) + pz*axesA(sp1,m1,6)
	tz = px*axesA(sp1,m1,7) + py*axesA(sp1,m1,8) + pz*axesA(sp1,m1,9)

	! Adjust position for different reference point on central molecule
	offx = tx - origin(1)
	offy = ty - origin(2)
	offz = tz - origin(3)

	! Check minimum / maximum separation between point and specified origin
	distsq = offx*offx + offy*offy + offz*offz
	if ((distsq.lt.mindistsq).or.(distsq.gt.maxdistsq)) then
	  getbin = .false.
	  return
	end if

	! Calculate integer position in grid (using non-offset coordinates)
	nx=NINT(offx/griddelta)
	ny=NINT(offy/griddelta)
	nz=NINT(offz/griddelta)

	if (MAX(ABS(nx),MAX(ABS(ny),ABS(nz))).gt.grid) then
	  getbin = .false.
	else
	  getbin = .true.
	end if

	end function getbin
