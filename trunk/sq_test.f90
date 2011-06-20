!	** sq **
!	Compute the static, neutron-weighted Faber-Ziman total structure factor S(Q)
!	Last changed 13 February 2009 - Partitioning of elements into 'types' implemented

	program sofq
	use dlprw; use utility; use parse
	implicit none
	include "/usr/local/mpich/include/mpif.h"
	!include "/usr/local/include/mpif.h"
	!include "/u/tristan/src/mpich-1.2.7p1/include/mpif.h"
	!include "/t/myri/mpich-gm/include/mpif.h"

	! Isotope definitions
	integer, parameter :: NISOTOPES = 9
	character*6 :: isonames(NISOTOPES)
	real*8 :: isoscatter(NISOTOPES), isofrac(NISOTOPES), isopops(NISOTOPES)
	integer :: isontypes(NISOTOPES), isonamelens(NISOTOPES)

	! DL_POLY (FF) to 'new type' to isotope mappings
	integer :: ntypes
	character*6, allocatable :: origtype(:), newtype(:), uniquetypes(:)
	integer, allocatable :: isotypes(:), typemap(:), namelens(:), uniqueiso(:)
	real*8, allocatable :: typepops(:), typefrac(:)

	! General variables
	real*8, parameter :: pi = 3.14159265358979
	character*80 :: hisfile,dlpoutfile,basename,resfile,namemap,headerfile
	character*20 :: temp
	integer :: i,j,k,baselen,bin,nframes,success,nargs,framestodo,nfftypes,alpha,beta
	integer :: n, m, o, nbins, found, kx, ky, kz, newkx, newky, newkz, nvec, newnvec, frameskip, sumfac, discard, framesdone
	integer :: nproc_mpi, id_mpi, err_mpi, mpistat(MPI_STATUS_SIZE)
	integer :: kpernode, kremain, kstart, kend, stdevmax, aoff, sp
	integer :: iargc
	real*8 :: binwidth,factor,weight,x1,x2,x3, totalweight, tx, ty, tz
	logical :: MASTER, SLAVE, writepartials = .FALSE., readmap = .FALSE., altheader = .FALSE.
	logical :: npt = .FALSE., calcstdev = .FALSE.
	integer, allocatable :: frameadded(:), kvectors(:,:), slaveadded(:), totaladded(:)
	real*8 :: kcut,magx,magy,magz,mag, numdensity, sd, sumxsq, avg
	real*8, allocatable :: framesq(:,:,:),sq(:),slavesq(:,:,:),partialsq(:,:,:)
	real*8, allocatable :: stdevdata(:,:,:,:), kmag(:)
	complex*16, allocatable :: rxxx(:,:),ryyy(:,:),rzzz(:,:),pdensity(:)
	character :: c

	! Scattering lengths for H,D,C,N,O,F,P,S,Cl
	isoscatter = (/ -3.706, 6.671, 6.646, 9.36, 5.803, 5.654, 5.13, 2.847, 9.577 /)
	isonames = (/ "H ", "D ", "C ", "N ", "O ", "F ", "P ", "S ", "Cl" /)
	isonamelens = (/ 1,1,1,1,1,1,1,1,2 /)

	binwidth=0.1   ! In Angstroms
	kcut = 5.0    ! Reciprocal space cutoff (box integers)
	frameskip = 1	! Take consecutive frames by default
	framestodo = -1	! Do all available frames by default
	discard = 0

	nargs = iargc()
	if (nargs.LT.2) then
	  write(0,*) "Usage : sq <DLP HISTORYfile> <DLP OUTPUTfile> ...options"
	  write(0,*) "            [-bin binwidth] [-frames nframes] [-kcut cutoff] [-discard nframes]"
	  write(0,*) "            [-discard interval] [-partials] [-readmap <file>] [-altheader <file>]"
	  write(0,*) "            [-stdev maxframes] [-npt]"
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
	      case ("-frames","-nframes"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
	      case ("-discard"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") discard
	      case ("-discard"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") frameskip
	      case ("-partials"); writepartials = .TRUE.
	      case ("-npt"); npt = .TRUE.
	      case ("-stdev"); calcstdev = .TRUE.; n = n + 1; call getarg(n,temp); read(temp,"(I6)") stdevmax
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
	  if (imcon.gt.2) stop "This image convention not properly supported - fix me!"

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
	  if (calcstdev) write(12,"(a,i6)") "Standard Deviations will be calculated per point, maxframes = ", stdevmax
	  if (outinfo(dlpoutfile,1).EQ.-1) then
	    write(12,*) "Problem with OUTPUT file."
	    call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    goto 999
	  end if

	  ! Read in forcefield names, the related internal types, and the corresponding isotope/element symbols
	  open(unit=15,file="lengths.dat",form="formatted",status="old")
	  read(15,"(I4)") nfftypes
	  allocate(origtype(nfftypes))
	  allocate(newtype(nfftypes))
	  allocate(isotypes(nfftypes))
	  allocate(typemap(natms))
	  isotypes(n) = 0
	  ! First column is original atom name, second is new type name, third is related element/isotope
	  do n=1,nfftypes
	    success = readline(15)
	    origtype(n) = arg(1)
	    newtype(n) = arg(2)
	    !read(15,"(A6,A6,)") origtype(n),newtype(n),temp
	    do m=1,NISOTOPES
	      if (arg(3).eq.isonames(m)) isotypes(n) = m
	    end do
	    if (isotypes(n).eq.0) then
	      write(0,*) "Unrecognised element/isotope in lengths.dat: ",temp
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
	      atmname(n) = newtype(found)
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
	    typemap(n) = found		! Store the type number for atom n
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
	if (SLAVE) allocate(typemap(natms))

	call MPI_BCast(typemap,natms,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(ntypes,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)

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
	    write(12,"(i2,2x,a6,1x,f7.5,3x,f7.4,2x,i4)") n,isonames(n),isofrac(n),isoscatter(n),isontypes(n)
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
	allocate(framesq(ntypes,ntypes,nbins))
	allocate(frameadded(nbins))
	allocate(pdensity(ntypes))

	! Allocate some arrays on the master
	if (MASTER) then
	  allocate(partialsq(ntypes,ntypes,nbins))
	  allocate(totaladded(nbins))
	  partialsq = 0.0
	  totaladded = 0
	  allocate(slavesq(ntypes,ntypes,nbins))
	  allocate(slaveadded(nbins))
	  if (calcstdev) then
	    if (MASTER) write(12,"(a,f10.4,a)") "Size of stdev array is :", ntypes*ntypes*nbins*stdevmax*4.0 / (1024.0*1024.0), "mb"
	    allocate(stdevdata(ntypes,ntypes,nbins,stdevmax))
	    stdevdata = 0.0
	  end if
	end if
	  
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
101	  success=readframe()
	  if (success.EQ.1) then   ! End of file encountered....
	    !open(unit=12,file=basename(1:baselen)//"out",form="formatted",status="old")
	    write(12,*) "End of history file found."
	    !close(12)
	    call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    goto 120
	  end if
	  if (success.EQ.-1) then  ! File error....
	    !open(unit=12,file=basename(1:baselen)//"out",form="formatted",status="old")
	    write(12,*) "History file ended prematurely..."
	    !close(12)
	    call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    goto 120
	  end if
	  nframes=nframes+1
	  numdensity = numdensity + natms / (cell(1)*cell(5)*cell(9))
	  !if (mod(nframes,100).EQ.0) then
	    !open(unit=12,file=basename(1:baselen)//"out",form="formatted",status="old")
	    write(6,*) nframes
	    !close(12)
	  !end if
	  if (nframes.LE.discard) goto 101	! Discard frames at beginning of trajectory if requested
	  if (mod(nframes,frameskip).ne.0) goto 101	! Frame skip interval

	  ! TEST - Reconstruct molecules straddling cell edges
	  ! First, calculate COM as out reference
	  call calc_com()
	  aoff = 0
	  do sp=1,nspecies
	    do m=1,s_nmols(sp)
	      do i=1,s_natoms(sp)
		call pbc(xpos(aoff+i),ypos(aoff+i),zpos(aoff+i),comx(sp,m),comy(sp,m),comz(sp,m),tx,ty,tz)
		xpos(aoff+i) = tx
		ypos(aoff+i) = ty
		zpos(aoff+i) = tz
	      end do
	      aoff = aoff + s_natoms(sp)
	    end do
	  end do
	  
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
	call MPI_BCast(kcut,1,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(npt,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_mpi)
	
	! Recalculate number of kvectors (if first frame or variable cell)
	if (npt.or.(framesdone.eq.0)) then

	  if (MASTER) call calc_rcell
	  call MPI_BCast(rcell,9,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	  newkx = (kcut / rcell(1)) + 1
	  newky = (kcut / rcell(5)) + 1
	  newkz = (kcut / rcell(9)) + 1
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
	        magx = i*rcell(1)
	        magy = j*rcell(5) 
	        magz = k*rcell(9)
	        mag = sqrt(magx**2 + magy**2 + magz**2)
	        if (mag.LT.kcut) newnvec = newnvec + 1
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
	          magx = i*rcell(1)
	          magy = j*rcell(5) 
	          magz = k*rcell(9)
	          mag = sqrt(magx**2 + magy**2 + magz**2)
	          if (mag.LT.kcut) then
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
	  write(13+id_mpi,"(a,i2,a,e14.6)") "Process ",id_mpi," thinks that kcut is ",kcut
	  write(13+id_mpi,"(a,i2,a,i7,a)") "Process ",id_mpi," thinks there are ",nproc_mpi," processes"
	  write(13+id_mpi,"(a,i2,a,i7,a)") "Process ",id_mpi," thinks there are ",nvec," vectors in total"
	  write(13+id_mpi,"(a,i2,a,i7,a,i7,a)") "Process ",id_mpi," calculating kvec=",kstart,",",kend,"..."
	end if

	! Each node constructs its own duplicate of the atomic k-vectors
	do n=1,natms
	  rxxx(0,n) = cmplx(1.0,0.0)
	  ryyy(0,n) = cmplx(1.0,0.0)
	  rzzz(0,n) = cmplx(1.0,0.0)
	  rxxx(1,n) = cmplx( cos(rcell(1)*xpos(n)) , sin(rcell(1)*xpos(n)) )
	  ryyy(1,n) = cmplx( cos(rcell(5)*ypos(n)) , sin(rcell(5)*ypos(n)) )
	  rzzz(1,n) = cmplx( cos(rcell(9)*zpos(n)) , sin(rcell(9)*zpos(n)) )
	  rxxx(-1,n) = conjg(rxxx(1,n))
	  ryyy(-1,n) = conjg(ryyy(1,n))
	  rzzz(-1,n) = conjg(rzzz(1,n))
	  do i=2,kx
	    rxxx(i,n) = rxxx(1,n) * rxxx(i-1,n)
	    rxxx(-i,n) = conjg(rxxx(i,n))
	  end do
	  do j=2,ky
	    ryyy(j,n) = ryyy(1,n) * ryyy(j-1,n)
	    ryyy(-j,n) = conjg(ryyy(j,n))
	  end do
	  do k=2,kz
	    rzzz(k,n) = rzzz(1,n) * rzzz(k-1,n)
	    rzzz(-k,n) = conjg(rzzz(k,n))
	  end do
	end do

	! Sum over reciprocal space vectors
	framesq = 0.0
	frameadded = 0
	do i=kstart,kend
	  pdensity = cmplx(0.0,0.0)
	  if (kvectors(i,1).eq.0) then
	    sumfac = 1
	  else
	    sumfac = 2
	  end if
	  do n=1,natms
	    ! Calculate p(r), sorting by atom type
	    alpha = typemap(n)
	    pdensity(alpha) = pdensity(alpha) + rxxx(kvectors(i,1),n) * ryyy(kvectors(i,2),n) * rzzz(kvectors(i,3),n)
	  end do
	  bin = kvectors(i,4)
	  frameadded(bin) = frameadded(bin) + sumfac
	  ! Calculate and store the partial structure factors S_ab(Q)
	  do alpha=1,ntypes
	    do beta=alpha,ntypes
	      framesq(alpha,beta,bin) = framesq(alpha,beta,bin) + real(sumfac) * (pdensity(alpha)*conjg(pdensity(beta)))
	      if (alpha.ne.beta) framesq(beta,alpha,bin) = framesq(beta,alpha,bin) + real(sumfac) * (pdensity(beta)*conjg(pdensity(alpha)))
	     !#if ((bin.eq.141).or.(bin.eq.155)) write(0,*) alpha,beta,bin,real(sumfac) * (pdensity(alpha)*conjg(pdensity(beta)))
	    end do
	  end do
	end do

	! Gather slave S(Q) data for this frame
	if (MASTER) then
	  ! Gather all partial S(Q) and frameadded data from slave processes into arrays on the master
	  write(0,*) "Master receiving slave data"
	  do i=1,nproc_mpi-1
	    call MPI_Recv(slavesq,ntypes*ntypes*nbins,MPI_REAL8,i,i+100,MPI_COMM_WORLD,mpistat,err_mpi)
	    ! Add this data into the local framesq array
	    do alpha=1,ntypes
	      do beta=1,ntypes
		do n=1,nbins
		  framesq(alpha,beta,n) = framesq(alpha,beta,n) + slavesq(alpha,beta,n)
		end do
	      end do
	    end do
	  end do
	  do i=1,nproc_mpi-1
	    call MPI_Recv(slaveadded,nbins,MPI_INTEGER,i,i+200,MPI_COMM_WORLD,mpistat,err_mpi)
	    do n=1,nbins
	      frameadded(n) = frameadded(n) + slaveadded(n)
	    end do
	  end do

	  ! Add to local partial sq
	  partialsq = partialsq + framesq
	  totaladded = totaladded + frameadded

	  ! Accumulate STDEV data
	  if (calcstdev) then
	    do alpha=1,ntypes
	      do beta=1,ntypes
		do n=1,nbins
		  stdevdata(alpha,beta,n,framesdone+1) = framesq(alpha,beta,n)
		end do
	      end do
	    end do
	  end if
	else
	  ! Slaves just send their data
	  write(0,*) "Slave sending data....", nbins
	  call MPI_Send(framesq,ntypes*ntypes*nbins,MPI_REAL8,0,id_mpi+100,MPI_COMM_WORLD,mpistat,err_mpi)
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

	  ! Normalisation of S_ab(Q)
	  do alpha=1,ntypes
	    do beta=1,ntypes
	      ! Atom type and bin populations
	      factor = dsqrt( typepops(alpha) * typepops(beta) )
	      do n=1,nbins
		if (totaladded(n).eq.0) then
		  partialsq(alpha,beta,n) = 0.0
		  if (calcstdev) stdevdata(alpha,beta,n,:) = 0
		else
		  partialsq(alpha,beta,n) = partialsq(alpha,beta,n) / (factor * totaladded(n))
		  if (calcstdev) stdevdata(alpha,beta,n,:) = stdevdata(alpha,beta,n,:) / (factor * totaladded(n))
		end if
	      end do
	      ! Normalise (shift) all partials to oscillate around 1.0 (will be uniformly subtracted later on)
	      do n=1,nbins
		weight = 1.0
		if (alpha.eq.beta) weight = 0.0
		partialsq(alpha,beta,n) = weight + partialsq(alpha,beta,n)
		if (calcstdev) stdevdata(alpha,beta,n,:) = weight + stdevdata(alpha,beta,n,:)
	      end do
	    end do
	  end do
	  write(0,*) "Finished normalise..."
	  
	  ! Sum into total S(Q)
	  totalweight = 0.0d0
	  sq = 0.0d0
	  do alpha=1,ntypes
	    do beta=1,ntypes

	      ! Weight against system density and uniformly subtract 1.0
	      do n=1,nbins
	        partialsq(alpha,beta,n) = (partialsq(alpha,beta,n) - 1.0) * 4.0 * pi * numdensity
	      end do

	      ! Write raw, unweighted partial S(Q) data if requested
	      if (writepartials) then
	        resfile=basename(1:baselen)//"rawsq"//uniquetypes(alpha)(1:namelens(alpha))//"-"//uniquetypes(beta)(1:namelens(beta))
	        OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	        do n=1,nbins
		  sd = 0.0
		  if (calcstdev) then
	            ! Calculate range values
	            avg = 0.0d0
	            do m=1,framesdone
	              avg = avg + stdevdata(alpha,beta,n,m)
	            end do
	            avg = avg / real(framesdone)
	            ! Calculate standard deviation
	            sumxsq = 0.0d0
	            do m=1,framesdone
	              sumxsq = sumxsq + (stdevdata(alpha,beta,n,m) - avg)**2
	            end do
	            sd = SQRT( sumxsq / real(framesdone) )
		  end if
		  c = " "
		  if (totaladded(n).eq.0) c = "#"
		  write(9,"(a,1x,F8.5,3x,E14.5,2x,e14.5,I8)") c, (n*binwidth)-binwidth*0.5, partialsq(alpha,beta,n), sd, totaladded(n)
	        end do
	        close(9)
	      end if

	      ! Weight accordingt to fractional populations
	      weight = typefrac(alpha)*typefrac(beta)
	      ! Additional weighting from partitioning of elements into 'sub-types'
	      weight = weight / dsqrt((typefrac(alpha)/isofrac(uniqueiso(alpha)))*(typefrac(beta)/isofrac(uniqueiso(beta))))
	      totalweight = totalweight + weight
	      write(6,"(a,f9.6,a,f9.6,a,f12.6)") "Partial "//uniquetypes(alpha)(1:namelens(alpha))//"-"//uniquetypes(beta)(1:namelens(beta))//" has fractional pops ",typefrac(alpha)," and ",typefrac(beta), "for a total atomic weighting of ", weight
	      ! Apply weighting to individual S_ab(Q)
	      do n=1,nbins
	        partialsq(alpha,beta,n) = weight * partialsq(alpha,beta,n)
	      end do

	      ! Weight by scattering length
	      weight = isoscatter(uniqueiso(alpha))*isoscatter(uniqueiso(beta))
	      write(6,"(a,f9.6,a,f9.6,a,f12.6)") "Partial "//uniquetypes(alpha)(1:namelens(alpha))//"-"//uniquetypes(beta)(1:namelens(beta))//" has scattering lengths ",isoscatter(uniqueiso(alpha))," and ", isoscatter(uniqueiso(beta)), "for a total neutron weight of ", weight
	      ! Apply weighting to individual S_ab(Q)
	      do n=1,nbins
	        partialsq(alpha,beta,n) = weight * partialsq(alpha,beta,n)
	      end do
	      ! Perform summation into total S(Q)
	      do n=1,nbins
	        sq(n) = sq(n) + partialsq(alpha,beta,n)
	      end do
	      ! Write partials if required
	      if (writepartials) then
	        resfile=basename(1:baselen)//"partsq"//uniquetypes(alpha)(1:namelens(alpha))//"-"//uniquetypes(beta)(1:namelens(beta))
	        OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	        do n=1,nbins
		  c = " "
		  if (totaladded(n).eq.0) c = "#"
		  write(9,"(a,1x,F8.5,3x,E14.5,2x,e14.5,I8)") c,(n*binwidth)-binwidth*0.5, partialsq(alpha,beta,n), sd, totaladded(n)
	        end do
	        close(9)
	      end if
	    end do
	  end do
	  write(6,*) "Finished sum"
	  write(6,*) "Total weight of summed partials:", totalweight

	  ! Calculate form constant (NOT USED)
	  factor = 0.0
	  do n=1,NISOTOPES
	    if (isopops(n).eq.0) cycle
	    factor = factor + isofrac(n) * isoscatter(n)
	  end do
	  factor = factor * factor

	  ! Write total S(Q)
	  resfile=basename(1:baselen)//"sq"
	  open(unit=9,file=resfile,form='formatted')
	  write(9,"(a,i8)") "# Un-normalised data follows. Total frames used was ",nframes
	  write(9,"(a,f10.5)") "# For information only, form constant (FZ) <b>**2 = ",factor

	  do n=1,nbins
	    c = " "
	    if (totaladded(n).eq.0) c = "#"
	    write(9,"(a,1x,F8.5,2x,E15.7,2x,I8)") c,(n-0.5)*binwidth,sq(n),totaladded(n)
	  end do
	  close(9)

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
		  write(9,"(f8.5,3x,E14.5,2x,I8)") (n*binwidth)-binwidth*0.5, sq(n)
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
