!	** sq **
!	Compute the static, neutron-weighted Faber-Ziman total structure factor S(Q)
!	Last changed 13 February 2009 - Partitioning of elements into 'types' implemented
!	Last changes 22 April 2010 - Added sanity check, calculating total scattering function without partitioning

	program sofq
	use dlprw; use utility; use parse
	implicit none
	include "mpif.h"

	! Isotope definitions
	integer, parameter :: NISOTOPES = 12
	character*6 :: isonames(NISOTOPES)
	real*8 :: isoscatter(NISOTOPES), isofrac(NISOTOPES), isopops(NISOTOPES)
	integer :: isontypes(NISOTOPES), isonamelens(NISOTOPES)

	! DL_POLY (FF) to 'new type' to isotope mappings
	integer :: ntypes, nexchangeable
	character*6, allocatable :: origtype(:), newtype(:), uniquetypes(:)
	integer, allocatable :: isotypes(:), typemap(:), namelens(:), uniqueiso(:), exchangelist(:)
	real*8, allocatable :: typepops(:), typefrac(:)
	logical, allocatable :: exchangeable(:)

	! General variables
	real*8, parameter :: pi = 3.14159265358979d0
	character*80 :: hisfile,dlpoutfile,basename,resfile,namemap,headerfile,lengthsfile
	character*20 :: temp
	integer :: i,j,k,baselen,bin,nframes,nargs,framestodo,nfftypes,alpha,beta
	integer :: n, m, o, nbins, found, kx, ky, kz, newkx, newky, newkz, nvec, newnvec, frameskip, sumfac, framestodiscard, framesdone
	integer :: nproc_mpi, id_mpi, err_mpi, mpistat(MPI_STATUS_SIZE)
	integer :: kpernode, kremain, kstart, kend
	integer :: iargc
	real*8 :: binwidth,factor,weight,x1,x2,x3, selfscatter, fac1, fac2
	logical :: MASTER, SLAVE, writepartials = .FALSE., readmap = .FALSE., altheader = .FALSE., npt = .FALSE.
	logical :: success
	integer, allocatable :: frameadded(:), kvectors(:,:), slaveadded(:), totaladded(:)
	real*8 :: kcut,kmin,magx,magy,magz,mag, numdensity, rposx, rposy, rposz
	real*8, allocatable :: framesq(:,:,:), sq(:), slavesq(:,:,:), partialsq(:,:,:), sanitysq(:), framesanitysq(:)
	real*8, allocatable :: stdev(:,:,:), m2n(:,:,:), kmag(:)
	complex*16, allocatable :: rxxx(:,:),ryyy(:,:),rzzz(:,:),pdensity(:)
	complex*16 :: sanitydensity, tempcomp
	character :: c

	! Scattering lengths for XX,H,D,C,N,O,F,P,S,Cl,Si
	isoscatter = (/ 0.0, -3.706, 6.671, 6.646, 9.36, 5.803, 5.654, 5.13, 2.847, 9.577, 4.1491, 7.718 /)
	isonames = (/ "X ", "H ", "D ", "C ", "N ", "O ", "F ", "P ", "S ", "Cl", "Si", "Cu" /)
	isonamelens = (/ 1,1,1,1,1,1,1,1,1,2,2,2 /)

	binwidth=0.1	! In Angstroms
	kcut = 5.0	! Reciprocal space cutoff (box multiples)
	kmin = 0.0	! Minimum cutoff
	frameskip = 1	! Take consecutive frames by default
	framestodo = -1	! Do all available frames by default
	framestodiscard = 0
        lengthsfile="lengths.dat"

	nargs = iargc()
	if (nargs.LT.2) then
	  write(0,*) "Usage : sq <HISTORYfile> <OUTPUTfile> ...options"
	  write(0,*) "            [-bin width]        Set x binwidth to use in S(Q) output"
	  write(0,*) "            [-frames n]         Set maximum number of frames to accumulate"
	  write(0,*) "            [-kcut cutoff]      Set maximum k magnitude to bin"
	  write(0,*) "            [-kmin cutoff]      Set minimum k magnitude to bin"
	  write(0,*) "            [-discard n]        Ignore n frames at the start of the trajectory"
	  write(0,*) "            [-skip n]           Skip n frames between each calculation"
	  write(0,*) "            [-partials]         Write partial S(Q) info"
	  write(0,*) "            [-readmap <file>]   Read alternative atom names map from <file>"
	  write(0,*) "            [-altheader <file>] Use specified DL_POLY history <file> for header"
	  write(0,*) "            [-npt]              Recalculate kvectors for each frame"
	  write(0,*) "            [-lengths <file>]   Use specified file instead of 'lengths.dat'"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	if (nargs.GE.3) then
	  n = 3
	  do
	    call getarg(n,temp)
	    select case (temp)
	      case ("-bin"); n = n + 1; call getarg(n,temp); read(temp,"(F20.10)") binwidth
	      case ("-kcut"); n = n + 1; call getarg(n,temp); read(temp,"(F20.10)") kcut
	      case ("-kmin"); n = n + 1; call getarg(n,temp); read(temp,"(F20.10)") kmin
	      case ("-frames","-nframes"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
	      case ("-discard"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodiscard
	      case ("-skip"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") frameskip
	      case ("-partials"); writepartials = .TRUE.
	      case ("-npt"); npt = .TRUE.
	      case ("-lengths"); n = n + 1; call getarg(n,lengthsfile)
	      case ("-readmap"); readmap = .TRUE.; n = n + 1; call getarg(n,namemap)
	      case ("-altheader"); altheader = .TRUE.; n = n + 1; call getarg(n,headerfile)
	      case ("-p4amslave")
		write(6,*) "Discarded option:",temp
	      case ("-p4wd","-p4pg","-execer_id","-master_host","-my_hostname","-my_nodenum","-my_numprocs","-total_numnodes","-master_port")
		write(6,*) "Discarded option:",temp
		n = n + 1; call getarg(n,temp); write(6,*) "   ...and argument:",temp
	      case ("-remote_info")
		write(6,*) "Discarded option: -remote_info"
		do
		  n = n + 1
		  if (n.gt.nargs) exit
		  call getarg(n,temp)
		  if (temp(1:1).eq."-") exit
		  write(6,*) "   ...and argument:",temp
		end do
	      case default
		write(6,*) "Unrecognised command line argument:", temp
		!stop
	    end select
	    n = n + 1
	    if (n.GT.nargs) exit
	  end do
	end if
	
	nbins = kcut / binwidth + 1

	! Initialise MPI and distribute data
	call MPI_Init(err_mpi)
	call MPI_Comm_rank(MPI_COMM_WORLD,id_mpi,err_mpi)
	call MPI_Comm_size(MPI_COMM_WORLD,nproc_mpi,err_mpi)

	if ((nproc_mpi.le.1).or.(id_mpi.eq.0)) then
	  MASTER = .true.; SLAVE = .false.
	else
	  MASTER = .false.; SLAVE = .true.
	end if
 
	if (MASTER) then
	  ! Open trajectory file
	  call openhis(hisfile,10)
	  if (readheader().EQ.-1) then
	    if (altheader) then
	      write(6,*) "Restarted trajectory:"
	      close(dlpun_his)
	      call openhis(headerfile,10)
	      if (readheader().EQ.-1) then
	        write(12,*) "Couldn't read header of alterhative history file."
	        call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	        goto 999
	      end if
	      close(dlpun_his)
	      call openhis(hisfile,10)
	    else
	      write(12,*) "Couldn't read header of history file."
	      goto 999
	    end if
	  end if

	  ! Check imcon
	  if (imcon.gt.3) stop "This image convention is not supported."

	  ! Ascertain length of basename....
	  baselen=-1
	  do n=80,1,-1
	    if (hisfile(n:n).EQ.".") THEN
	      baselen=n
	      goto 2
	    endif
	  end do
2         if (baselen.EQ.-1) THEN
	    basename="sq."
	    baselen=6
	  else
	    basename=hisfile(1:baselen)
	  endif

	  open(unit=12,file=basename(1:baselen)//"out",form="formatted",status="replace")
	  write(12,*) "There are ",nproc_mpi," MPI processes."
	  write(12,"(A,F6.3,A)") "Using binwidth of ",binwidth," Angstroms"
	  write(12,"(A,I5,A)") "There will be ",nbins," histogram bins."
	  if (outinfo(dlpoutfile,1).EQ.-1) then
	    write(12,*) "Problem with OUTPUT file."
	    call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    goto 999
	  end if

	  ! Read in forcefield names, the related internal types, and the corresponding isotope/element symbols
	  open(unit=15,file=lengthsfile,form="formatted",status="old")
	  read(15,"(I4)") nfftypes
	  allocate(origtype(nfftypes))
	  allocate(newtype(nfftypes))
	  allocate(isotypes(nfftypes))
	  allocate(exchangeable(nfftypes))
	  allocate(typemap(natms), exchangelist(natms))
	  isotypes = 0
	  ! First column is original atom name, second is new type name, third is related element/isotope, fourth is exchangeable flag
	  do n=1,nfftypes
	    success = readline(15)
	    origtype(n) = arg(1)
	    newtype(n) = arg(2)
	    write(0,"(A6,A6,a6)") origtype(n),newtype(n),arg(3)
	    do m=1,NISOTOPES
	      if (arg(3).eq.isonames(m)) isotypes(n) = m
	    end do
	    if (isotypes(n).eq.0) then
	      write(0,*) "Unrecognised element/isotope in lengths.dat: ",arg(3)
	      call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	      goto 999
	    end if
	    ! Check that this type only 'uses' a single isotope
	    do m=1,n-1
	      if ((newtype(n).eq.newtype(m)).and.(isotypes(n).ne.isotypes(m))) then
		write(0,"(a,a,a,a,a)") "Error - mapping original type ",origtype(n)," to new type ", newtype(n)," uses a different isotope than one assigned earlier"
		call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
		goto 999
	      end if
	    end do
	    ! Is this an exchangeable hydrogen?
	    if ((arg(3)(1:1).eq."x").or.(arg(3)(1:1).eq."X")) then
	      ! Check that this is H or D
	      if ((arg(3)(1:1).ne."H").and.(arg(3)(1:1).ne."D")) then
		write(0,"(a,a)") "Error - exchangeable flag set for non-H/D : ",arg(3)
		call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
		goto 999
	      end if
	      ! Set the flag
	      exchangeable(n) = .TRUE.
	    else
	      exchangeable(n) = .FALSE.
	    end if
	  end do
	  close(15)

	  ! If specified, read in a different set of atom names from a specified file (-readmap)
	  ! Allows for the changing of atom names or easier splitting of FF types to 1 or more types
	  if (readmap) then
	    write(12,"(A,A)") "Reading alternative list of atom names from: ",namemap
	    open(unit=15,file=namemap,form="formatted",status="old",err=50)
	    do n=1,natms
	      read(15,"(A6)",err=50) atmname(n)
	    end do
	    goto 55
50	    write(12,"(A,I5,A,I5)") "Error reading from namefile for atom ",n," of ",natms
55	    close(15)
	    stop "Failed."
	  end if

	  ! Now that the atomnames->partitions->isotopes map has been defined, we convert HISfile atomnames into partition names
	  ! At the same time, construct a unique list of new type names
	  allocate(uniquetypes(nfftypes))
	  allocate(uniqueiso(nfftypes))
	  uniquetypes = ""
	  uniqueiso = 0
	  ntypes = 0
	  do n=1,natms
	    ! First, convert atom name from forcefield type in history file to element type specified in 'lengths.dat'
	    found = -1
	    do m=1,nfftypes
	      if (origtype(m).eq.atmname(n)) found = m
	    end do
	    if (found.eq.-1) then
	      write(0,"(a,a,a)") "Atomtype ",atmname(n)," not listed in lengths.dat."
	      call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	      goto 999
	    else
	      write(12,"(a,i6,a,a,a,i2,a,a)") "Atom ",n,", typename ",atmname(n)," matches entry ",found, " in lengths.dat : ",origtype(found)
	      atmname(n) = newtype(found)
	      write(12,"(a,a)") "...mapped to new name ",atmname(n)
	    end if
	    ! Now compare the name of the atom with those stored in 'typenames' to see if it's 'new'
	    found = -1
	    do m=1,ntypes
	      if (uniquetypes(m).eq.atmname(n)) found = m
	    end do
	    ! Add this type to the list if it has not been encountered before
	    if (found.EQ.-1) then
	      ntypes = ntypes + 1
	      uniquetypes(ntypes) = atmname(n)
	      found = ntypes
	    end if
	    ! Store the unique type ID for the atom
	    typemap(n) = found
	  end do

	  ! Assign isotypes to the unique types list
	  isontypes = 0
	  do n=1,ntypes
	    uniqueiso(n) = 0
	    do m=1,nfftypes
	      if (uniquetypes(n).eq.newtype(m)) uniqueiso(n) = isotypes(m)
	    end do
	    if (isotypes(n).eq.0) then
	      write(12,*) "Unrecognised element/isotope in unique list - this shouldn't have happened! : ",uniquetypes(n)
	      call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	      goto 999
	    end if
	    isontypes(uniqueiso(n)) = isontypes(uniqueiso(n)) + 1
	  end do
	  write(12,"(A,I2,A)") "Atoms are to be partitioned into ",ntypes," unique groups:"

	  ! Send out proceed flag to slaves
	  call MPI_BCast(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	else
	  ! Slaves just wait for flag
	  call MPI_BCast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  if (i.ne.1) goto 999
	end if
	  
	call MPI_BCast(natms,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(nfftypes,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	if (SLAVE) then
	  allocate(typemap(natms))
	  allocate(uniqueiso(nfftypes))
	end if

	call MPI_BCast(typemap,natms,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(ntypes,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(uniqueiso,nfftypes,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)

	if (MASTER) then
	  allocate(typepops(ntypes))
	  allocate(typefrac(ntypes))
	  allocate(namelens(ntypes))
	  ! Calculate the populations of the new atom types *and* fractional populations of related isotopes
	  typepops = 0
	  typefrac = 0
	  isofrac = 0
	  isopops = 0
	  do n=1,natms
	    typepops(typemap(n)) = typepops(typemap(n)) + 1.0
	    isopops(uniqueiso(typemap(n))) = isopops(uniqueiso(typemap(n))) + 1
	  end do
	  isofrac = isopops / natms
	  typefrac = typepops / natms

	  write(12,*) "   Type      N      Frac   El/Iso   S.Length"
	  write(12,"(i2,2x,a6,1x,f7.0,3x,f7.4,2x,a6,2x,f7.4)") (n,uniquetypes(n),typepops(n),typepops(n)/natms,isonames(uniqueiso(n)),isoscatter(uniqueiso(n)),n=1,ntypes)

	  write(12,*) ""
	  write(12,*) "Element/isotope fractional populations:"
	  write(12,*) "  Name     Frac    S.Length   NTypes"
	  do n=1,NISOTOPES
	    if (isopops(n).eq.0) cycle
	    write(12,"(i2,2x,a6,1x,f8.5,3x,f7.4,2x,i4)") n,isonames(n),isofrac(n),isoscatter(n),isontypes(n)
	  end do

	  ! Set the name lengths of each of the typenames
	  do n=1,ntypes
	    do m=6,1,-1
	      if (uniquetypes(n)(m:m).EQ." ") namelens(n) = m-1
	    end do
	  end do
	
	  deallocate(newtype)
	  deallocate(origtype)
	end if

	! Slaves must manually allocate their own position arrays, and needs the bin value from the master
	call MPI_BCast(nbins,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	if (SLAVE) then
	  allocate(xpos(natms))
	  allocate(ypos(natms))
	  allocate(zpos(natms))
	end if

	allocate(sq(nbins))
	allocate(framesanitysq(nbins))
	allocate(framesq(ntypes,ntypes,nbins))
	allocate(frameadded(nbins))
	allocate(pdensity(ntypes))

	! Allocate some arrays on the master
	if (MASTER) then
	  allocate(partialsq(ntypes,ntypes,nbins))
	  allocate(sanitysq(nbins))
	  allocate(totaladded(nbins))
	  partialsq = 0.0
	  sanitysq = 0.0
	  totaladded = 0
	  allocate(slavesq(ntypes,ntypes,nbins))
	  allocate(slaveadded(nbins))
	  allocate(stdev(ntypes,ntypes,nbins),m2n(ntypes,ntypes,nbins))
	  stdev = 0.0
	  m2n = 0.0
	end if

	call MPI_BCast(kcut,1,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(kmin,1,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(npt,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_mpi)
	  
	! MAIN LOOP BEGINS

	nframes=0
	framesdone = 0
	kx = 0
	ky = 0
	kz = 0
	nvec = 0
	numdensity = 0.0

100	if (MASTER) then
	  ! The master will read in the configuration and send the coords out to the slaves
101	  n = readframe()
	  if (n.EQ.1) then   ! End of file encountered....
	    !open(unit=12,file=basename(1:baselen)//"out",form="formatted",status="old")
	    write(12,*) "End of history file found."
	    !close(12)
	    call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    goto 120
	  end if
	  if (n.EQ.-1) then  ! File error....
	    !open(unit=12,file=basename(1:baselen)//"out",form="formatted",status="old")
	    write(12,*) "History file ended prematurely..."
	    !close(12)
	    call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    goto 120
	  end if
	  nframes=nframes+1
	  numdensity = numdensity + natms / volume(cell)
	  !if (mod(nframes,100).EQ.0) then
	    !open(unit=12,file=basename(1:baselen)//"out",form="formatted",status="old")
	    write(6,*) nframes
	    !close(12)
	  !end if
	  if (nframes.LE.framestodiscard) goto 101	! Discard frames at beginning of trajectory if requested
	  if (mod(nframes,frameskip).ne.0) goto 101	! Frame skip interval
	  ! Send out proceed 'flag' to slaves
	  call MPI_BCast(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	else
	  ! Slaves just wait for flag
	  call MPI_BCast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  if (i.EQ.0) goto 120
	end if

	call MPI_BCast(xpos,natms,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(ypos,natms,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(zpos,natms,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	
	! Recalculate number of kvectors (if first frame or variable cell)
	if (npt.or.(framesdone.eq.0)) then

	  if (MASTER) call calc_rcell()
	  call MPI_BCast(rcell,9,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	  newkx = (kcut / dsqrt(rcell(1)*rcell(1) + rcell(2)*rcell(2) + rcell(3)*rcell(3))) + 1
	  newky = (kcut / dsqrt(rcell(4)*rcell(4) + rcell(5)*rcell(5) + rcell(6)*rcell(6))) + 1
	  newkz = (kcut / dsqrt(rcell(7)*rcell(7) + rcell(8)*rcell(8) + rcell(9)*rcell(9))) + 1
	  if (MASTER) write(6,"(a,f10.6,a,3i4)") "This frame: rcell(1) = ",rcell(1)," with kxyz = ", newkx, newky, newkz
	  ! Reallocate if newkx is greater than kx
	  if ((newkx.gt.kx).or.(newky.gt.ky).or.(newkz.gt.kz)) then
	    kx = newkx
	    ky = newky
	    kz = newkz
	    if (framesdone.gt.0) then
	      deallocate(rxxx)
	      deallocate(ryyy)
	      deallocate(rzzz)
	    end if
	    allocate(rxxx(-kx:kx,natms))
	    allocate(ryyy(-ky:ky,natms))
	    allocate(rzzz(-kz:kz,natms))
	    if (MASTER) write(0,"(a)") "--> Reallocated working arrays."
	  end if

	  ! Calculate the new number of kvectors that will result from the choice of cutoff.
	  newnvec = 0
	  do i=0,kx
	    do j=-ky,ky
	      do k=-kz,kz
		magx = i*rcell(1) + j*rcell(2) + k*rcell(3)
		magy = i*rcell(4) + j*rcell(5) + k*rcell(6)
		magz = i*rcell(7) + j*rcell(8) + k*rcell(9)
	        mag = sqrt(magx**2 + magy**2 + magz**2)
	        if ((mag.le.kcut).and.(mag.ge.kmin)) newnvec = newnvec + 1
	      end do
	    end do
	  end do
	  if (MASTER) write(12,"(a,i8,a)") "This frame needs space for ",newnvec," k-vectors."
	  if (newnvec.gt.nvec) then
	    nvec = newnvec
	    if (framesdone.gt.0) deallocate(kvectors)
	    allocate(kvectors(nvec,4))
	    nvec = 0
	    do i=0,kx
	      do j=-ky,ky
	        do k=-kz,kz
		  magx = i*rcell(1) + j*rcell(2) + k*rcell(3)
		  magy = i*rcell(4) + j*rcell(5) + k*rcell(6)
		  magz = i*rcell(7) + j*rcell(8) + k*rcell(9)
	          mag = sqrt(magx**2 + magy**2 + magz**2)
	          if ((mag.le.kcut).and.(mag.ge.kmin)) then
 		    nvec = nvec + 1
 		    kvectors(nvec,1) = i
		    kvectors(nvec,2) = j
		    kvectors(nvec,3) = k
		    kvectors(nvec,4) = mag / binwidth + 1
	          end if
	        end do
	      end do
	    end do
	  end if
	  if (MASTER) write(12,"(a)") "--> Reallocated kvectors() array."
	  ! Work out new kvector ranges for nodes
	  kpernode = nvec / nproc_mpi
	  kremain = mod(nvec,kpernode)
	  kstart = id_mpi*kpernode + 1
	  kend = kstart + kpernode - 1
	  if (id_mpi+1.eq.nproc_mpi) kend = kend + kremain
	  write(13+id_mpi,"(a,i2,a,e14.6,a,e14.6)") "Process ",id_mpi," thinks that kcut is ",kcut, "and kmin is",kmin
	  write(13+id_mpi,"(a,i2,a,i7,a)") "Process ",id_mpi," thinks there are ",nproc_mpi," processes"
	  write(13+id_mpi,"(a,i2,a,i7,a)") "Process ",id_mpi," thinks there are ",nvec," vectors in total"
	  write(13+id_mpi,"(a,i2,a,i7,a,i7,a)") "Process ",id_mpi," calculating kvec=",kstart,",",kend,"..."
	end if

	! Each node constructs its own duplicate of the atomic k-vectors
	do n=1,natms
	  rxxx(0,n) = cmplx(1.0,0.0)
	  ryyy(0,n) = cmplx(1.0,0.0)
	  rzzz(0,n) = cmplx(1.0,0.0)
	  rposx = xpos(n)*rcell(1) + ypos(n)*rcell(4) + zpos(n)*rcell(7)
	  rposy = xpos(n)*rcell(2) + ypos(n)*rcell(5) + zpos(n)*rcell(8)
	  rposz = xpos(n)*rcell(3) + ypos(n)*rcell(6) + zpos(n)*rcell(9)
	  rxxx(1,n) = cmplx( cos(rposx) , sin(rposx) )
	  ryyy(1,n) = cmplx( cos(rposy) , sin(rposy) )
	  rzzz(1,n) = cmplx( cos(rposz) , sin(rposz) )
	  rxxx(-1,n) = dconjg(rxxx(1,n))
	  ryyy(-1,n) = dconjg(ryyy(1,n))
	  rzzz(-1,n) = dconjg(rzzz(1,n))
	  do i=2,kx
	    rxxx(i,n) = rxxx(1,n) * rxxx(i-1,n)
	    rxxx(-i,n) = dconjg(rxxx(i,n))
	  end do
	  do j=2,ky
	    ryyy(j,n) = ryyy(1,n) * ryyy(j-1,n)
	    ryyy(-j,n) = dconjg(ryyy(j,n))
	  end do
	  do k=2,kz
	    rzzz(k,n) = rzzz(1,n) * rzzz(k-1,n)
	    rzzz(-k,n) = dconjg(rzzz(k,n))
	  end do
	end do

	! Sum over reciprocal space vectors
	framesq = 0.0
	framesanitysq = 0.0
	frameadded = 0
	do i=kstart,kend
	  if (MASTER.and.(mod(i-kstart,1000).eq.0)) write(12,"('... ',f4.1,' completed')") real(i-kstart)/real(kend-kstart) * 100.0
	  pdensity = cmplx(0.0,0.0)
	  sanitydensity = cmplx(0.0,0.0)
	  if (kvectors(i,1).eq.0) then
	    sumfac = 1
	  else
	    sumfac = 2
	  end if
	  do n=1,natms
	    ! Calculate p(r), sorting by atom type
	    alpha = typemap(n)
	    tempcomp = rxxx(kvectors(i,1),n) * ryyy(kvectors(i,2),n) * rzzz(kvectors(i,3),n) 
	    pdensity(alpha) = pdensity(alpha) + tempcomp
	    sanitydensity = sanitydensity + tempcomp * isoscatter(uniqueiso(alpha))
	  end do
	  bin = kvectors(i,4)
	  frameadded(bin) = frameadded(bin) + sumfac
	  ! Calculate and store the partial structure factors S_ab(Q)
	  do alpha=1,ntypes
	    do beta=alpha,ntypes
	      framesq(alpha,beta,bin) = framesq(alpha,beta,bin) + real(sumfac) * (pdensity(alpha)*dconjg(pdensity(beta)))
	      if (alpha.ne.beta) framesq(beta,alpha,bin) = framesq(beta,alpha,bin) + real(sumfac) * (pdensity(beta)*dconjg(pdensity(alpha)))
	    end do
	  end do
	  framesanitysq(bin) = framesanitysq(bin) + real(sumfac) * (sanitydensity*dconjg(sanitydensity))
	end do

	! Gather slave S(Q) data for this frame
	if (MASTER) then
	  ! Store local sanity S(Q) data
	  sanitysq = sanitysq + framesanitysq
	  ! Gather all partial S(Q) and frameadded data from slave processes into arrays on the master
	  !write(0,*) "Master receiving slave data"
	  do i=1,nproc_mpi-1
	    ! Receive partial S(Q) data
	    call MPI_Recv(slavesq,ntypes*ntypes*nbins,MPI_REAL8,i,i+100,MPI_COMM_WORLD,mpistat,err_mpi)
	    ! Add this data into the local framesq array
	    do alpha=1,ntypes
	      do beta=1,ntypes
		do n=1,nbins
		  framesq(alpha,beta,n) = framesq(alpha,beta,n) + slavesq(alpha,beta,n)
		end do
	      end do
	    end do
	    ! Receive sanity S(Q) data from slaves
	    call MPI_Recv(framesanitysq,nbins,MPI_REAL8,i,i+100,MPI_COMM_WORLD,mpistat,err_mpi)
	    sanitysq = sanitysq + framesanitysq
	  end do
	  do i=1,nproc_mpi-1
	    call MPI_Recv(slaveadded,nbins,MPI_INTEGER,i,i+200,MPI_COMM_WORLD,mpistat,err_mpi)
	    do n=1,nbins
	      frameadded(n) = frameadded(n) + slaveadded(n)
	    end do
	  end do

	  ! Accumulate STDEV data
	  if (framesdone.gt.0) then
	    fac1 = 1.0d0 / real(framesdone+1)
	    fac2 = 1.0d0 / real(framesdone)
	    m2n = m2n + (framesq - fac1*(framesq+partialsq))*(framesq - fac2*partialsq)
	  end if

	  ! Add to local partial sq
	  partialsq = partialsq + framesq
	  totaladded = totaladded + frameadded

	else
	  ! Slaves just send their data
	  !write(0,*) "Slave sending data....", nbins
	  call MPI_Send(framesq,ntypes*ntypes*nbins,MPI_REAL8,0,id_mpi+100,MPI_COMM_WORLD,mpistat,err_mpi)
	  call MPI_Send(framesanitysq,nbins,MPI_REAL8,0,id_mpi+100,MPI_COMM_WORLD,mpistat,err_mpi)
	  call MPI_Send(frameadded,nbins,MPI_INTEGER,0,id_mpi+200,MPI_COMM_WORLD,mpistat,err_mpi)
	end if

	framesdone = framesdone + 1
	if (MASTER) then
	  !open(unit=12,file=basename(1:baselen)//"out",form="formatted",status="old")
	  write(12,*) "Frames completed : ",framesdone
	  !close(12)
	end if

	! Next frame
	if (framesdone.eq.framestodo) goto 120
	goto 100

120	if (MASTER) then

	  write(12,*) "Finished."
	  write(12,*) ""

	  ! Get average number density
	  numdensity = numdensity / nframes
	  write(0,*) "Average atomic number density = ",numdensity

	  ! Sum into total S(Q)
	  sq = 0.0d0
	  selfscatter = 0.0
	  do alpha=1,ntypes
	    do beta=1,ntypes

	      ! Calculate weighting factor for partial (does not include fractional population since this is accounted for)
	      factor = isoscatter(uniqueiso(alpha)) * isoscatter(uniqueiso(beta))

	      ! Increment self-scattering correction
	      weight = 0.0
	      if (alpha.eq.beta) weight = typefrac(alpha) * factor
	      selfscatter = selfscatter + weight

	      ! Normalisation of S_ab(Q) w.r.t. number of frames and (total) number of atoms
	      do n=1,nbins
		if (totaladded(n).eq.0) then
		  partialsq(alpha,beta,n) = 0.0
		else
		  partialsq(alpha,beta,n) = partialsq(alpha,beta,n) / totaladded(n) / natms
		  m2n(alpha,beta,n) = dsqrt(m2n(alpha,beta,n) / totaladded(n) / natms)
		end if
		sq(n) = sq(n) + partialsq(alpha,beta,n) * factor
	      end do
	      write(0,*) "Finished normalise..."

	      ! Write partials if required
	      if (writepartials) then
	        resfile=basename(1:baselen)//"partsq"//uniquetypes(alpha)(1:namelens(alpha))//"-"//uniquetypes(beta)(1:namelens(beta))
	        OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
		write(9,"('# b_a * b_b = ',f10.5)") factor
	        write(9,"('# Self-scattering level is ',f10.5)") weight
		write(9,"(a68)") "#  Q           S_ab(Q)        STDEV       S_ab(Q)*bibj   NAdded   S_ab(Q)-Self"
	        do n=1,nbins
		  write(9,"(F8.5,3x,E14.5,2x,e14.5,1x,e14.5,1x,I10,1x,e14.5)") (n-0.5)*binwidth, partialsq(alpha,beta,n), m2n(alpha,beta,n), partialsq(alpha,beta,n)*factor, totaladded(n), partialsq(alpha,beta,n)-weight
	        end do
	        close(9)
	      end if
	    end do
	  end do
	  write(6,*) "Finished sum"

	  ! Remove self-scattering from S(Q) and convert between units of fm (www.ncnr.nist.gov) and 
	  ! 10**-12 cm units (ISIS) used for the scattering lengths
	  sq = (sq - selfscatter) / 100.0

	  ! Normalise sanity S(Q)
	  do n=1,nbins
	    sanitysq(n) = sanitysq(n) / totaladded(n) / natms
	  end do
	  !sanitysq = sanitysq * 4.0 * pi * numdensity / 100.0
	  sanitysq = (sanitysq - selfscatter) / 100.0

	  ! Calculate form constant (NOT USED)
	  factor = 0.0
	  do n=1,NISOTOPES
	    if (isopops(n).eq.0) cycle
	    factor = factor + isofrac(n) * isoscatter(n) * isoscatter(n)
	  end do
	  factor = factor * factor

	  ! Write total S(Q) and sanity S(Q)
	  resfile=basename(1:baselen)//"sq"
	  open(unit=9,file=resfile,form='formatted')
	  write(9,"(a,i8)") "# F(Q) data follows (ds/dS - selfscatter). Total frames used was ",nframes
	  write(9,"(a,f10.5)") "# Self scattering (sum (c * b**2)) = ", selfscatter
	  write(9,"(a,f10.5)") "# For information only, form constant (FZ) <b>**2 = ",factor
	  write(9,"(a,f10.5)") "# For information only, 4*PI*rho = ", 4.0 * pi * numdensity
	  resfile=basename(1:baselen)//"sanitysq"
	  open(unit=8,file=resfile,form='formatted')
	  write(8,"(a,i8)") "# F(Q) data follows (ds/dS - selfscatter). Total frames used was ",nframes
	  write(8,"(a,f10.5)") "# Self scattering (sum (c * b**2)) = ", selfscatter
	  write(8,"(a,f10.5)") "# For information only, form constant (FZ) <b>**2 = ",factor
	  write(8,"(a,f10.5)") "# For information only, 4*PI*rho = ", 4.0 * pi * numdensity

	  do n=1,nbins
	    c = " "
	    if (totaladded(n).eq.0) c = "#"
	    write(9,"(a,1x,F8.5,2x,E15.7,2x,I10)") c,(n-0.5)*binwidth,sq(n),totaladded(n)
	    write(8,"(a,1x,F8.5,2x,E15.7,2x,I10)") c,(n-0.5)*binwidth,sanitysq(n),totaladded(n)
	  end do
	  close(9)
	  close(8)

	  ! Write elemental (isotopic) partials if required
	  if (writepartials) then
	    write(6,*) "Writing partials..."
	    do i=1,NISOTOPES
	      if (isontypes(i).eq.0) cycle
	      do j=i,NISOTOPES
	        if (isontypes(j).eq.0) cycle
	        ! We will re-use the S(Q) array
	        sq = 0.0
	        do alpha=1,ntypes
	          do beta=1,ntypes
		    if (((uniqueiso(alpha).eq.i).and.(uniqueiso(beta).eq.j)).or.((uniqueiso(alpha).eq.j).and.(uniqueiso(beta).eq.i))) then
		      do n=1,nbins
			sq(n) = sq(n) + partialsq(alpha,beta,n)
		      end do
		    end if
		  end do
		end do
		! Write data
		resfile=basename(1:baselen)//"elemsq"//isonames(i)(1:isonamelens(i))//"-"//isonames(j)(1:isonamelens(j))
		OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
		do n=1,nbins
		  write(9,"(f8.5,3x,E14.5)") (n*binwidth)-binwidth*0.5, sq(n)
		end do
		close(9)
	      end do
	    end do
	  end if

	  close(12)
	  write(6,*) "Master finished!"
	else
	  write(6,*) "Slave ",id_mpi," finished."
	end if

999	call MPI_Finalize(err_mpi)
	end program sofq
