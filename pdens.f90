!	** pdens.f90 **
!	Program to calculate the 3-dimensional distribution of different species about
!	one species' centres.
!	Updated version - rotates all molecules about target molecule first, then performs minimum
!	image before binning

	program calcpdens
	use dlprw; use utility
	implicit none
	integer, parameter :: MAXSITES = 20, MAXSPECIES = 5
	real*8, parameter :: radcon = 57.29577951d0
	real*8, allocatable :: pdens(:,:,:,:)
	real*8, allocatable :: avggeom(:,:)		! Average species coordinates
	real*8, allocatable :: pdensintra(:,:,:)	! Intramolecular species distributions
	integer, allocatable :: nfound(:), ncaught(:)	! Total and 'binned' mols
	integer, allocatable :: molflags(:)		! Per-frame map of sp1 mols to include in averaging
	integer, allocatable :: spexp(:)		! Expected species numbers
	integer :: atomSites(MAXSPECIES,0:MAXSITES)	! Atoms for surrounding species centres (if any)
	integer :: grid = 20, pgrid = 50		! Grid in each direction for 3ddist and pdens
	real*8 :: delta = 0.5, pdelta = 0.15		! Default grid spacings
	logical :: molmap = .FALSE.			! Whether we're using a file of mol flags
	logical :: nopdens = .FALSE., nointra = .FALSE.	! Flags to prohibit various calculations
	integer :: nmapframes, baselen,nargs,nframesused
	character*80 :: hisfile,outfile,basename,resfile,temp,flagfile,altheaderfile
	logical :: altheader = .FALSE.
	integer :: success
	integer :: n, nframes, m, o, p, i
	real*8 :: px, py, pz, mindistsq, maxdistsq, distsq, tx, ty, tz
	real*8 :: orientangle(20,3), orientdelta(20,3), angdelta, ax
	integer :: centresp,sp1,m1,sp2,m2,startf,endf,npdenscentres,nintracentres,nx,ny,nz
	logical :: orientcheck(20,3), failed, getbin
	integer :: iargc

	write(0,*) "*** pdens"

	nargs = iargc()
	if (nargs.lt.3) then
	  write(*,"(a)") "Usage: pdens <HIStory file> <OUTput file> <centresp> [-options]"
	  write(*,"(a)") "        [-axis sp x1 x2 y1 y2]     Atoms to use for axis calculation in species sp"
	  write(*,"(a)") "        [-atom sp i]               Add atom i in species sp as a specific atom to use as other point"
	  write(*,"(a)") "        [-delta spacing]           Spacing between grid points (default = 0.5 Angstrom)"
	  write(*,"(a)") "        [-end frameno]             Trajectory frame to end calculations on (default = last)"
	  write(*,"(a)") "        [-grid npoints]            Grid points in each direction for prob. densities (default = 20)"
	  write(*,"(a)") "        [-header file]             Use specified file to get header"
	  write(*,"(a)") "        [-maxdist r]               Maximum separation allowed between molecules (default = 10000.0)"
	  write(*,"(a)") "        [-mindist r]               Minimum separation allowed between molecules (default = 0.0)"
	  write(*,"(a)") "        [-molmap <file>]           Formatted file (I5,1X,n(I2,1x) of frame,mols to use of sp1"
	  write(*,"(a)") "        [-nointer]                 Prohibit calculation of the intermolecular probability densities"
	  write(*,"(a)") "        [-nointra]                 Prohibit calculation of the intramolecular probability density and average molecule"
	  write(*,"(a)") "        [-orient sp axis angle delta]  Specify a restriction in the angular orientation of the second &
		& molecule (use negative delta for symmetry about 90deg)"
	  write(*,"(a)") "        [-otheraxis sp x1 x2 y1 y2] Alternative axis definition, for surrounding &
		& molecule position and orientation"
	  write(*,"(a)") "        [-pdelta spacing]          Spacing between pgrid points (default = 0.15 Angstrom)"
	  write(*,"(a)") "        [-pgrid npoints]           Grid points in each direction for species densities (default = 50)"
	  write(*,"(a)") "        [-start frameno]           Trajectory frame to start calculations (default = 1)"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
	  call getarg(3,temp); read (temp,"(i4)") centresp
	end if
	n = 10
	write(0,"(A,A)") "History file : ",hisfile
	!call openhis(hisfile,n)
	write(0,"(A,A)") " Output file : ",outfile
	if (outinfo(outfile,1).eq.-1) goto 798
	write(0,"(A,i3)") "Central species is = ",centresp

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
	mindistsq = 0.0
	maxdistsq = 1.0e9
	centresp = 1
	atomSites = 0
	startf = 1
	endf = 0
	npdenscentres = 0
	nintracentres = 0
	orientcheck = .FALSE.
	orientangle = 0.0
	orientdelta = 0.0

	n = 3
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-axis") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp1
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp1,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp1,2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp1,3)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp1,4)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp1," calculated from: X=",axesAatoms(sp1,1),"->", &
	        & axesAatoms(sp1,2),", Y=0.5X->0.5(r(",axesAatoms(sp1,3),")->r(",axesAatoms(sp1,4),"))"
	      axesAdefined(sp1) = .true.
	    case ("-atom")
	      n = n + 1; call getarg(n,temp); read(temp,"(I4)") sp1
	      if (axesBAtoms(sp1,1).ne.0) stop "Definition of sites and otheraxis cannot be used together."
	      atomSites(sp1,0) = atomSites(sp1,0) + 1
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") atomSites(sp1, atomSites(sp1,0))
	      write(0,"(A,I3,a,i2)") "Atom ", atomSites(sp1, atomSites(sp1,0)), " added as other site for species ", sp1
	    case ("-grid")
	      n = n + 1; call getarg(n,temp); read(temp,"(I4)") grid
	      write(0,"(A,I4)") "Grid points in each XYZ = ",grid
	    case ("-pgrid")
	      n = n + 1; call getarg(n,temp); read(temp,"(I4)") pgrid
	      write(0,"(A,I4)") "PGrid points in each XYZ = ",pgrid
	    case ("-delta")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") delta
	      write(0,"(A,F6.3)") "Grid spacing = ",delta
	    case ("-pdelta")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") pdelta
	      write(0,"(A,F6.3)") "PGrid spacing = ",pdelta
	    case ("-start")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") startf
	      write(0,"(A,I5)") "Starting frame = ",startf
	    case ("-end")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") endf
	      write(0,"(A,I5)") "End frame = ",endf
	    case ("-molmap")
	      n = n + 1; call getarg(n,flagfile)
	      write(0,"(A,A)") "Using flags for species 1 molecules from : ",flagfile
	      open(unit=20,file=flagfile,form="formatted",status="old")
	      molmap = .TRUE.
	    case ("-nopdens")
	      nopdens = .TRUE.
	      write(0,"(A,I4)") "Intermolecular probability distributions will not be computed."
	      write(15,"(A,I4)") "Intermolecular probability distributions will not be computed."
	    case ("-nointra")
	      nointra = .TRUE.
	      write(0,"(A,I4)") "Intramolecular probability distributions will not be computed."
	      write(15,"(A,I4)") "Intramolecular probability distributions will not be computed."
	    case ("-mindist")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") mindistsq
	      write(0,"(A,f8.4)") "Minimum separation between molecules = ",mindistsq
	      write(15,"(A,f8.4)") "Minimum separation between molecules = ",mindistsq
	      mindistsq = mindistsq * mindistsq
	    case ("-maxdist")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") maxdistsq
	      write(0,"(A,f8.4)") "Maximum separation between molecules = ",maxdistsq
	      write(15,"(A,f8.4)") "Maximum separation between molecules = ",maxdistsq
	      maxdistsq = maxdistsq * maxdistsq
	    case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
              altheader = .TRUE.
	    case ("-orient")
	      n = n + 1; call getarg(n,temp); read(temp,"(i10)") sp1
	      n = n + 1; call getarg(n,temp); read(temp,"(i10)") i
	      n = n + 1; call getarg(n,temp); read(temp,"(f20.14)") orientangle(sp1,i)
	      n = n + 1; call getarg(n,temp); read(temp,"(f20.14)") orientdelta(sp1,i)
	      write(0,"(A,i2,a,i1,a,f10.4,a,f10.4,a)") "Only molecules of sp ",sp1," with axis ", i, " angle delta of ", &
		& orientdelta(sp1,i), " degrees about ", orientangle(sp1,i), " will be binned"
	      orientcheck(sp1,i) = .TRUE.
	    case ("-otheraxis") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp1
	      if (atomSites(sp1,0).ne.0) stop "Definition of sites and otheraxis cannot be used together."
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesBatoms(sp1,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesBatoms(sp1,2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesBatoms(sp1,3)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesBatoms(sp1,4)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Alternative axes for species ",sp1," calculated from: X=",axesBatoms(sp1,1),"->", &
	        & axesBatoms(sp1,2),", Y=0.5X->0.5(r(",axesBatoms(sp1,3),")->r(",axesBatoms(sp1,4),"))"
	      axesBdefined(sp1) = .true.
	    case default
	      write(0,*) "Unrecognised command line option : ",temp
	      stop
	  end select
	end do
	
	! Print out a summary of the control variables to the output file.
	write(15,"(A,A)") "Input file: ",hisfile
	write(15,"(A,I3)") "Molecular species in file : ",nspecies
	do sp1=1,nspecies
	  ! Check that the necessary molecules have had their axes defined
	  if (axesAdefined(sp1)) then
	    write(15,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp1," calculated from: X=",axesAatoms(sp1,1),"->", &
	      & axesAatoms(sp1,2),", Y=0.5X->0.5(r(",axesAatoms(sp1,3),")->r(",axesAatoms(sp1,4),"))"
	  else if (sp1.eq.centresp) then
	    stop "A proper set of axes must be defined for the central species."
	  else if (orientcheck(sp1,1).or.orientcheck(sp1,2).or.orientcheck(sp1,3)) then
	    stop "Axes must be defined on any species whose orientation is to be considered."
	  end if
	end do
	write(15,"(A,I4)") "Grid (intermolecular) points in each XYZ = ",grid
	write(15,"(A,I4)") "PGrid (intramolecular) points in each XYZ = ",pgrid
	write(15,"(A,F6.3)") "Grid spacing = ",delta
	write(15,"(A,F6.3)") "PGrid spacing = ",pdelta
	write(15,"(A,I4)") "Central species = ",centresp
	write(15,"(A,I5,A,I5,A)") "Frame range = ",startf," to ",endf," (0=last)"
	if (nointra.and.nopdens) stop "Nothing to do (both -nointra and -nopdens given)."

	! Probability density arrays
	allocate(pdens(nspecies,-grid:grid,-grid:grid,-grid:grid))
	allocate(ncaught(nspecies))
	allocate(nfound(nspecies))
	allocate(pdensintra(-pgrid:pgrid,-pgrid:pgrid,-pgrid:pgrid))
	allocate(spexp(nspecies))
	! Average species etc.
	allocate(avggeom(maxatoms,3))
	if (molmap) allocate(molflags(s_nmols(centresp)))

	! Clear the arrays before we start...
	write(0,*) "Clearing initial arrays..."
	pdens = 0.0
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
	    write(0,"(I6,2x,'(',I6,',',I10,',',I10')')") nframes,nmapframes,npdenscentres,nintracentres
	  else
	    write(0,"(i6,2x,'(',I10,',',I10,')')") nframes,npdenscentres,nintracentres
	  end if
	end if
	if (nframes.LT.startf) goto 102

	if ((molmap).and.(m1.ne.nframes)) goto 102
	if (molmap) nmapframes = nmapframes + 1

	nframesused = nframesused + 1

	! Generate all molecular axes
	call genaxes()

	!
	! Calculate 3D distributions
	!
	if (nopdens) goto 115
	sp1 = centresp
	do m1=1,s_nmols(sp1)      ! Loop over all molecules of species 1...
	  ! If we're using a mapfile, decide whether to include this molecule
	  if (molmap) then
	    if (molflags(m1).eq.0) cycle 
	  end if
	  npdenscentres = npdenscentres + 1	! Normalisation counter
	  do sp2=1,nspecies
	    p = s_start(sp2) - 1
	    do m2=1,s_nmols(sp2)

	      ! Don't consider central molecule with itself
	      if ((sp1.eq.sp2).and.(m1.eq.m2)) cycle

	      ! Check orientation if it has been specified
	      failed = .false.
	      do n=1,3
	        if (orientcheck(sp2,n)) then
	          ! Take dot product of X-axis on sp2 with that of the central molecule, and calculate the angle delta
		  if (axesBdefined(sp2)) then
		    ax = axesB(sp2,m2,(n-1)*3+1)*axesA(sp1,m1,(n-1)*3+1) + axesB(sp2,m2,(n-1)*3+2)*axesA(sp1,m1,(n-1)*3+2) + axesB(sp2,m2,(n-1)*3+3)*axesA(sp1,m1,(n-1)*3+3)
		  else
		    ax = axesA(sp2,m2,(n-1)*3+1)*axesA(sp1,m1,(n-1)*3+1) + axesA(sp2,m2,(n-1)*3+2)*axesA(sp1,m1,(n-1)*3+2) + axesA(sp2,m2,(n-1)*3+3)*axesA(sp1,m1,(n-1)*3+3)
		  end if
	          angdelta = acos(ax) * radcon
	          ! If delta is negative, map delta onto +/-90
	          if ((orientdelta(sp2,n).lt.0.0).and.(angdelta.gt.90.0)) angdelta = angdelta - 180.0
	          if (dabs(angdelta-orientangle(sp2,n)).gt.dabs(orientdelta(sp2,n))) failed = .true.
		end if
		if (failed) exit
	      end do
	      if (failed) cycle

	      ! Pass position of axesAorigin, or axesBorogin, or individual atomSites, depending on what has been defined
	      if (axesBdefined(sp2)) then
		if (getbin(sp1,m1,mindistsq,maxdistsq,axesBorigin(sp2,m2,1),axesBorigin(sp2,m2,2),axesBorigin(sp2,m2,3),nx,ny,nz,grid,delta)) then
		  pdens(sp2,nx,ny,nz) = pdens(sp2,nx,ny,nz) + 1
		  ncaught(sp2) = ncaught(sp2) + 1
		end if
		nfound(sp2) = nfound(sp2) + 1
	      else if (atomSites(sp2,0).gt.0) then
		if (getbin(sp1,m1,mindistsq,maxdistsq,axesAorigin(sp2,m2,1),axesAorigin(sp2,m2,2),axesAorigin(sp2,m2,3),nx,ny,nz,grid,delta)) then
		  pdens(sp2,nx,ny,nz) = pdens(sp2,nx,ny,nz) + 1
		  ncaught(sp2) = ncaught(sp2) + 1
		end if
		nfound(sp2) = nfound(sp2) + 1
	      else
		do n=1,atomSites(sp2,0)
		  i = p+atomSites(sp2,n)
		  if (getbin(sp1,m1,mindistsq,maxdistsq,xpos(i),ypos(i),zpos(i),nx,ny,nz,grid,delta)) then
		    pdens(sp2,nx,ny,nz) = pdens(sp2,nx,ny,nz) + 1
		    ncaught(sp2) = ncaught(sp2) + 1
		  end if
		  nfound(sp2) = nfound(sp2) + 1
		end do
	      end if

	      p=p+s_natoms(sp2)
	    end do
	  end do
	end do

	! Calculate average molecule and intramolecular distribution
115	if (nointra) goto 116
	p=s_start(centresp)
	do m1=1,s_nmols(centresp)
	  ! If we're using a mapfile, decide whether to include this molecule
	  if (molmap) then
	    if (molflags(m1).eq.0) then
	      p = p+s_natoms(centresp)
	      cycle
	    end if
	    nintracentres = nintracentres + 1
	  else
	    nintracentres = nintracentres + 1
	  end if
	  do n=1,s_natoms(centresp)
	    ! PBC the atom
	    call pbc(xpos(p+n-1),ypos(p+n-1),zpos(p+n-1),axesAorigin(centresp,m1,1), &
	  	& axesAorigin(centresp,m1,2),axesAorigin(centresp,m1,3),tx,ty,tz)
	    ! 'Zero' its position with respect to the axis centre...
	    tx=tx-axesAorigin(centresp,m1,1)
	    ty=ty-axesAorigin(centresp,m1,2)
	    tz=tz-axesAorigin(centresp,m1,3)
	    ! Transform the coordinate into the local coordinate system...
	    px=tx*axesA(centresp,m1,1) + ty*axesA(centresp,m1,2) + tz*axesA(centresp,m1,3)
	    py=tx*axesA(centresp,m1,4) + ty*axesA(centresp,m1,5) + tz*axesA(centresp,m1,6)
	    pz=tx*axesA(centresp,m1,7) + ty*axesA(centresp,m1,8) + tz*axesA(centresp,m1,9)
	    ! Accumulate the position....
	    avggeom(n,1) = avggeom(n,1)+px
	    avggeom(n,2) = avggeom(n,2)+py
	    avggeom(n,3) = avggeom(n,3)+pz
	    ! Intramolecular distribution.....
	    nx=NINT(px/pdelta)
	    ny=NINT(py/pdelta)
	    nz=NINT(pz/pdelta)
	    if (MAX(ABS(nx),MAX(ABS(ny),ABS(nz))).lt.pgrid) then
	      pdensintra(nx,ny,nz) = pdensintra(nx,ny,nz) + 1
	    else
	      write(0,*) "Large distance : nx,ny,nz = ",nx,ny,nz
	    end if
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

	write(0,"(a,i1,a,i7,a,i7,a)") "Intermolecular calculation consideres ",npdenscentres," molecules"
	write(0,"(a,i1,a,i7,a,i7,a)") "Intramolecular calculation consideres ",nintracentres," molecules"

	! Set normalisation over frames and central species molecules
	do sp1=1,nspecies

	  ! Calculate expected species numbers
	  spexp(sp1) = s_nmols(centresp) * s_nmols(sp1) * nframesused * max(1,atomSites(sp1,0))
	  if (centresp.eq.sp1) spexp(sp1) = spexp(sp1) - s_nmols(sp1) * nframesused * max(1,atomSites(sp2,0))

	  ! Species density about central species
	  pdens(sp1,:,:,:) = pdens(sp1,:,:,:)/npdenscentres/(delta**3)
	  !do nx=-grid,grid
	  !  do ny=-grid,grid
	  !    do nz=-grid,grid
	  !      pdens(sp1,nx,ny,nz)=pdens(sp1,nx,ny,nz)/n3dcentres/(delta**3)
	  !    end do
	  !  end do
	  !end do

	end do
 
	! Intramolecular distribution
        pdensintra(:,:,:) = pdensintra(:,:,:)/nintracentres/(pdelta**3)
	!do nx=-pgrid,pgrid
	!  do ny=-pgrid,pgrid
	!    do nz=-pgrid,pgrid
        !      pdensintra(sp1,nx,ny,nz)=pdensintra(sp1,nx,ny,nz)/(nframesused*s_nmols(sp1))/(pdelta**3)
	!    end do
	!  end do
	!end do

	! Average molecule
	avggeom(:,:) = avggeom(:,:)/nintracentres
	!do n=1,s_natoms(sp1)
	!  do m=1,3
	!    avggeom(sp1,n,m)=avggeom(sp1,n,m)/(nframesused*s_nmols(sp1))
	!  end do
	!end do

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
	do sp1=1,nspecies
	  write(0,"(A,I1,A,A)") " Species ",sp1,", ",s_name(sp1)
	  write(0,"(A12,I12,A,I9,A)") "Expected : ",spexp(sp1)," over ",nframesused," frames."
	  write(0,"(A12,I12)") "Found : ",nfound(sp1)
	  write(0,"(A12,I12,A)") "Caught : ",ncaught(sp1)," (in grid)"
	  if (sp1.eq.centresp) write(0,"(A12,I9)") "Selected (inter) : ",npdenscentres
	  if (sp1.eq.centresp) write(0,"(A12,I9)") "Selected (intra) : ",nintracentres
	  write(15,"(A,I1,A,A)") " Species ",sp1,", ",s_name(sp1)
	  write(15,"(A12,I12,A,I9,A)") "Expected : ",spexp(sp1)," over ",nframesused," frames."
	  write(15,"(A12,I12)") "Found : ",nfound(sp1)
	  write(15,"(A12,I12,A)") "Caught : ",ncaught(sp1)," (in grid)"
	  if (sp1.eq.centresp) write(15,"(A12,I9)") "Selected (inter) : ",npdenscentres
	  if (sp1.eq.centresp) write(15,"(A12,I9)") "Selected (intra) : ",nintracentres
	end do
	write(0,"(3(A,I3),A,F4.1,A)") "Grid = ",(grid*2+1),"x",(grid*2+1),"x", &
	  & (grid*2+1)," ( x ",delta," Angstrom )"
	write(15,"(3(A,I3),A,F4.1,A)") "Grid = ",(grid*2+1),"x",(grid*2+1),"x", &
	  & (grid*2+1)," ( x ",delta," Angstrom )"
	write(0,"(A)") "-- Writing output files..."

	do sp1=1,nspecies

	  if (.not.nopdens) then
	    ! -- Probability density of this molecule about the central species
	    resfile=basename(1:baselen)//CHAR(48+centresp)//CHAR(48+sp1)//".pdens"
	    open(unit=9,file=resfile,form="formatted")
	    write(9,*) 2*grid+1, 2*grid+1, 2*grid+1
	    write(9,"(9f8.4)") delta,0.0,0.0,0.0,delta,0.0,0.0,0.0,delta
	    write(9,"(3f10.4)") -grid*delta,-grid*delta,-grid*delta
	    write(9,*) "zyx"
	    do n=-grid,grid
	      do m=-grid,grid
	        do o=-grid,grid
	    	write(9,"(f12.8)") pdens(sp1,n,m,o)
	        end do
	      end do
	    end do
	    close(9)
	  end if

	end do

	if (.not.nointra) then
	  ! -- Average molecule
	  resfile=basename(1:baselen)//"avg"//char(48+centresp)//".xyz"
	  open(unit=9,file=resfile,form="formatted")
	  ! Write out pdb files of the average molecules
900	  FORMAT (a2,5x,3(F9.5,2x))
	  write(9,*) s_natoms(centresp)
	  write(9,*) "Average: ", s_name(centresp)
	  do n=1,s_natoms(centresp)
	    write(9,900) atmname(s_start(centresp)-1+n)(1:2),(avggeom(n,m),m=1,3)
	  end do
	  close(9)

	  ! -- Intramolecular probability distribution
	  resfile=basename(1:baselen)//"intra"//char(48+centresp)//".pdens"
	  open(unit=9,file=resfile,form="formatted")
	  write(9,*) 2*pgrid+1, 2*pgrid+1, 2*pgrid+1
	  write(9,"(9f8.4)") pdelta,0.0,0.0,0.0,pdelta,0.0,0.0,0.0,pdelta
	  write(9,"(3f10.4)") -pgrid*pdelta,-pgrid*pdelta,-pgrid*pdelta
	  write(9,*) "zyx"
	  do n=-pgrid,pgrid
	    do m=-pgrid,pgrid
	      do o=-pgrid,pgrid
	        write(9,"(F12.8)") pdensintra(n,m,o)
	      end do
	    end do
	  end do
	  close(9)
	end if

	write(0,*) "Finished!"
	write(15,"(A)") "Finished!"
999	close(10)
	close(13)

	end program calcpdens

	logical function getbin(sp1,m1,mindistsq,maxdistsq,x,y,z,nx,ny,nz,grid,griddelta)
	use utility
	implicit none
	integer, intent(in) :: sp1, m1, grid
	integer, intent(out) :: nx, ny, nz
	real*8, intent(in) :: x, y, z, mindistsq, maxdistsq
	real*8 :: px, py, pz, tx, ty, tz, distsq, griddelta

	! Check minimum / maximum separation
	distsq = x*x + y*y + z*z
	if ((distsq.lt.mindistsq).or.(distsq.gt.maxdistsq)) then
	  getbin = .false.
	  return
	end if

	! Get mim of supplied coordinate with the origin of the central species
	call pbc(x,y,z,axesAorigin(sp1,m1,1),axesAorigin(sp1,m1,2),axesAorigin(sp1,m1,3),tx,ty,tz)

	! Get vector between minimuim image coordinates and central molecule axis origin
	px = tx - axesAorigin(sp1,m1,1)
	py = ty - axesAorigin(sp1,m1,2)
	pz = tz - axesAorigin(sp1,m1,3)
	tx = px*axesA(sp1,m1,1) + py*axesA(sp1,m1,2) + pz*axesA(sp1,m1,3)
	ty = px*axesA(sp1,m1,4) + py*axesA(sp1,m1,5) + pz*axesA(sp1,m1,6)
	tz = px*axesA(sp1,m1,7) + py*axesA(sp1,m1,8) + pz*axesA(sp1,m1,9)

	! Calculate integer position in grid
	nx=NINT(tx/griddelta)
	ny=NINT(ty/griddelta)
	nz=NINT(tz/griddelta)

	if (MAX(ABS(nx),MAX(ABS(ny),ABS(nz))).gt.grid) then
	  getbin = .false.
	else
	  getbin = .true.
	end if

	end function getbin
