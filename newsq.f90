!	** sq **
!	Compute the static, neutron-weighted Faber-Ziman total structure factor S(Q)
!	Changed 16 March 2009 - Added alternative binning, fixed some MPI bugs.
!	Changed 13 February 2009 - Partitioning of elements into 'types' implemented

	program sofq
	use dlprw; use utility; use parse
	implicit none
	!include "/usr/local/mpich-1.2.7p1/include/mpif.h"
	!include "/usr/local/include/mpif.h"
	include "mpif.h"
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
	character*80 :: hisfile,dlpoutfile,basename,resfile,namemap
	character*20 :: temp
	integer :: i,j,k,baselen,bin,nframes,success,nargs,framestodo,nfftypes,alpha,beta
	integer :: n, m, o, nbins, found, kx, nvec, frameskip, sumfac, discard, framesdone
	integer :: nproc_mpi, id_mpi, err_mpi, mpistat(MPI_STATUS_SIZE)
	integer :: kpernode, kremain, kstart, kend
	integer :: iargc, altbins = 0
	real*8 :: binwidth,boxvolume,factor,weight,x1,x2,x3
	logical :: MASTER, SLAVE, writepartials = .FALSE., readmap = .FALSE.
	integer, allocatable :: numadded(:), kvectors(:,:), slaveadded(:)
	real*8 :: kcut,magx,magy,magz,mag, numdensity
	real*8, allocatable :: partialsq(:,:,:),sq(:),slavesq(:,:,:)
	complex*16, allocatable :: rxxx(:,:),ryyy(:,:),rzzz(:,:),pdensity(:)

	! Scattering lengths for H,D,C,N,O,F,P,S,Cl
	isoscatter = (/ -3.706, 6.671, 6.646, 9.36, 5.803, 5.654, 5.13, 2.847, 9.577 /)
	isonames = (/ "H ", "D ", "C ", "N ", "O ", "F ", "P ", "S ", "Cl" /)
	isonamelens = (/ 1,1,1,1,1,1,1,1,2 /)

	binwidth=0.1   ! In Angstroms
	kcut = 5.0    ! Reciprocal space cutoff (box integers)
	frameskip = 1	! Take consecutive frames by default
	discard = 0

	nargs = iargc()
	if (nargs.LT.2) then
	  write(0,*) "Usage : sq <DLP HISTORYfile> <DLP OUTPUTfile> ...options"
	  write(0,*) "                                [-bin binwidth] [-frames nframes] [-kcut cutoff] [-discard nframes]"
	  write(0,*) "                                [-skip interval] [-partials] [-readmap <file>] [-altbins]"
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
	      case ("-frames"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
	      case ("-discard"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") discard
	      case ("-skip"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") frameskip
	      case ("-partials"); writepartials = .TRUE.
	      case ("-altbins"); altbins = 1
	      case ("-readmap"); readmap = .TRUE.; n = n + 1; call getarg(n,namemap)
	      case default
		write(0,*) "Unrecognised command line argument:", temp
	    end select
	    n = n + 1
	    if (n.GT.nargs) exit
	  end do
	end if
	

	! Initialise MPI and distribute data
	call MPI_Init(err_mpi)
	call MPI_Comm_rank(MPI_COMM_WORLD,id_mpi,err_mpi)
	call MPI_Comm_size(MPI_COMM_WORLD,nproc_mpi,err_mpi)

	if (id_mpi.EQ.0) then
	  MASTER = .true.; SLAVE = .false.
	else
	  MASTER = .false.; SLAVE = .true.
	end if
 
	if (MASTER) then
	  ! Open trajectory file
	  call openhis(hisfile,10)

	  ! CHeck imcon
	  if (imcon.gt.1) stop "This image convention not properly supported - fix me!"

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
	  if (altbins) then
	    write(12,"(A)") "Using alternative quantized binning."
	  else
	    write(12,"(A,F6.3,A)") "Using binwidth of ",binwidth," Angstroms"
	  end if
	  if (outinfo(dlpoutfile,1).EQ.-1) then
	    write(12,*) "Problem with OUTPUT file."
	    goto 999
	  end if
	  ! Now, read in the history header so that we have cell() and also the atom names
	  if (readheader().EQ.-1) then
	    write(12,*) "Couldn't read header of history file."
	    goto 999
	  end if

	  ! Read in forcefield names, the related internal types, and the corresponding isotope/element symbols
	  open(unit=15,file="lengths.dat",form="formatted",status="old")
	  read(15,"(I4)") nfftypes
	  allocate(origtype(nfftypes))
	  allocate(newtype(nfftypes))
	  allocate(isotypes(nfftypes))
	  allocate(typemap(natms))
	  ! First column is original atom name, second is new type name, third is related element/isotope
	  do n=1,nfftypes
	    call readline(15)
	    origtype(n) = arg(1)
	    newtype(n) = arg(2)
	    !read(15,"(A6,A6,)") origtype(n),newtype(n),temp
	    isotypes(n) = 0
	    do m=1,NISOTOPES
	      if (arg(3).eq.isonames(m)) isotypes(n) = m
	    end do
	    if (isotypes(n).eq.0) then
	      write(12,*) "Unrecognised element/isotope in lengths.dat: ",temp
	      stop
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
	      stop
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
	      stop
	    end if
	    isontypes(uniqueiso(n)) = isontypes(uniqueiso(n)) + 1
	  end do
	  write(12,"(A,I2,A)") "Atoms are to be partitioned into ",ntypes," unique groups:"
	end if
	  
	call MPI_BCast(natms,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(framestodo,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(altbins,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
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

	! Calculate the reciprocal cell
	if (MASTER) call calc_rcell
	call MPI_BCast(rcell,9,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	call MPI_BCast(kcut,1,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	kx = (kcut / rcell(1)) + 1
	allocate(rxxx(-kx:kx,natms))
	allocate(ryyy(-kx:kx,natms))
	allocate(rzzz(-kx:kx,natms))
	if (altbins) then
	  nbins = 3*kx*kx + 1
	else
	  nbins = kcut / binwidth + 1
	end if
	if (MASTER) write(12,"(A,I5,A)") "There will be ",nbins," bins."
	write(0,"(a,i3,a,i5,a)") "Process ",id_mpi," thinks there are ", nbins, "bins"

	! Slaves must manually allocate their own position arrays
	if (SLAVE) then
	  allocate(xpos(natms))
	  allocate(ypos(natms))
	  allocate(zpos(natms))
	end if

	allocate(sq(0:nbins))
	allocate(partialsq(ntypes,ntypes,0:nbins))
	allocate(numadded(0:nbins))
	allocate(pdensity(ntypes))

	partialsq = 0.0
	numadded = 0

	! Calculate the number of kvectors that will result from the choice of cutoff.
	! Then, calculate all the vectors and send out to processes
	nvec = 0
	do i=0,kx
	  do j=-kx,kx
	    do k=-kx,kx
	      magx = i*rcell(1)
	      magy = j*rcell(5) 
	      magz = k*rcell(9)
	      mag = sqrt(magx**2 + magy**2 + magz**2)
	      if (mag.LT.kcut) nvec = nvec + 1
	    end do
	  end do
	end do
	if (MASTER) write(12,*) "There are ",nvec," k-vectors per frame."
	write(0,*) nvec
	allocate(kvectors(nvec,4))
	nvec = 0
	do i=0,kx
	  do j=-kx,kx
	    do k=-kx,kx
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
		!write(0,"(5i6)") nvec,i,j,k,i*i+j*j+k*k
	      end if
	    end do
	  end do
	end do

	! MAIN LOOP BEGINS

	kpernode = nvec / nproc_mpi
	kremain = mod(nvec,kpernode)
	nframes=0
	framesdone = 0

	!if (MASTER) close(12)

	! Work out the kvector range for each node
	kstart = id_mpi*kpernode + 1
	kend = kstart + kpernode - 1
	if (id_mpi+1.eq.nproc_mpi) kend = kend + kremain
	write(0,*) "Process ",id_mpi," calculating kvec=",kstart,",",kend,"..."

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
	    write(0,*) nframes
	    !close(12)
	  !end if
	  if (nframes.LE.discard) goto 101	! Discard frames at beginning of trajectory if requested
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
	
	! Each node constructs its own duplicate of the atomic k-vectors
	do n=1,natms
	  rxxx(0,n) = cmplx(1.0,0.0)
	  ryyy(0,n) = cmplx(1.0,0.0)
	  rzzz(0,n) = cmplx(1.0,0.0)
	  rxxx(1,n) = cmplx( cos(rcell(1)*xpos(n)) , sin(rcell(1)*xpos(n)) )
	  ryyy(1,n) = cmplx( cos(rcell(1)*ypos(n)) , sin(rcell(1)*ypos(n)) )
	  rzzz(1,n) = cmplx( cos(rcell(1)*zpos(n)) , sin(rcell(1)*zpos(n)) )
	  rxxx(-1,n) = conjg(rxxx(1,n))
	  ryyy(-1,n) = conjg(ryyy(1,n))
	  rzzz(-1,n) = conjg(rzzz(1,n))
	  do i=2,kx
	    rxxx(i,n) = rxxx(1,n) * rxxx(i-1,n)
	    ryyy(i,n) = ryyy(1,n) * ryyy(i-1,n)
	    rzzz(i,n) = rzzz(1,n) * rzzz(i-1,n)
	    rxxx(-i,n) = conjg(rxxx(i,n))
	    ryyy(-i,n) = conjg(ryyy(i,n))
	    rzzz(-i,n) = conjg(rzzz(i,n))
	  end do
	end do

	! Sum over reciprocal space vectors
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
	  if (altbins) then
	    bin = kvectors(i,1)*kvectors(i,1) + kvectors(i,2)*kvectors(i,2) + kvectors(i,3)*kvectors(i,3)
	  else
	    bin = kvectors(i,4)
	  endif
	  numadded(bin) = numadded(bin) + sumfac
	  ! Calculate and store the partial structure factors S_ab(Q)
	  do alpha=1,ntypes
	    do beta=alpha,ntypes
	      partialsq(alpha,beta,bin) = partialsq(alpha,beta,bin) + real(sumfac) * (pdensity(alpha)*conjg(pdensity(beta)))
	      if (alpha.ne.beta) partialsq(beta,alpha,bin) = partialsq(beta,alpha,bin) + real(sumfac) * (pdensity(beta)*conjg(pdensity(alpha)))
	    end do
	  end do
	end do

	framesdone = framesdone + 1
	write(0,*) "Process ",id_mpi," finished frame ",framesdone
	if (MASTER) then
	  !open(unit=12,file=basename(1:baselen)//"out",form="formatted",status="old")
	  write(12,*) "Frames completed : ",framesdone
	  !close(12)
	end if

	! Next frame
	if (framesdone.eq.framestodo) goto 120
	goto 100

120	if (MASTER) then
	  !open(unit=12,file=basename(1:baselen)//"out",form="formatted",status="old")
	  write(12,*) "Finished."
	  write(12,*) ""
	  ! Gather all partial S(Q) and numadded data from slave processes into arrays on the master
	  allocate(slavesq(ntypes,ntypes,nbins))
	  do i=1,nproc_mpi-1
	    write(0,*) "Receiving data from slave ",i, nbins
	    call MPI_Recv(slavesq,ntypes*ntypes*nbins,MPI_REAL8,i,i+100,MPI_COMM_WORLD,mpistat,err_mpi)
	    ! Add this data into the local partialsq array
	    do alpha=1,ntypes
	      do beta=1,ntypes
		do n=1,nbins
		  partialsq(alpha,beta,n) = partialsq(alpha,beta,n) + slavesq(alpha,beta,n)
		end do
	      end do
	    end do
	  end do
	  deallocate(slavesq)
	  allocate(slaveadded(nbins))
	  do i=1,nproc_mpi-1
	    call MPI_Recv(slaveadded,nbins,MPI_INTEGER,i,i+200,MPI_COMM_WORLD,mpistat,err_mpi)
	    do n=1,nbins
	      numadded(n) = numadded(n) + slaveadded(n)
	    end do
	  end do
	  deallocate(slaveadded)
	else
	  ! Slaves just send their data
	  write(0,*) "Slave sending data...", nbins
	  call MPI_Send(partialsq,ntypes*ntypes*nbins,MPI_REAL8,0,id_mpi+100,MPI_COMM_WORLD,mpistat,err_mpi)
	  call MPI_Send(numadded,nbins,MPI_INTEGER,0,id_mpi+200,MPI_COMM_WORLD,mpistat,err_mpi)
	end if

	if (MASTER) then

	  ! Get average number density
	  numdensity = numdensity / nframes
	  write(0,*) "Average atomic number density = ",numdensity

	  ! Normalisation of S_ab(Q)
	  do alpha=1,ntypes
	    do beta=1,ntypes
	      ! Atom type and bin populations
	      factor = dsqrt( typepops(alpha) * typepops(beta) )
	      do n=1,nbins
		if (numadded(n).eq.0) then
		  partialsq(alpha,beta,n) = 0.0
		else
		  partialsq(alpha,beta,n) = partialsq(alpha,beta,n) / (factor * numadded(n))
		end if
	      end do
	      ! Normalise (shift) all partials to oscillate around 1.0 (will be uniformly subtracted later on)
	      do n=1,nbins
		weight = 1.0
		if (alpha.eq.beta) weight = 0.0
		partialsq(alpha,beta,n) = weight + partialsq(alpha,beta,n)
	      end do
	    end do
	  end do
	  write(0,*) "Finished normalisation."
	  
	  ! Sum into total S(Q)
	  sq = 0.0d0
	  do alpha=1,ntypes
	    do beta=1,ntypes
	      ! Weight factor attributable to scattering length and fractional population
	      weight = typefrac(alpha)*isoscatter(uniqueiso(alpha))*typefrac(beta)*isoscatter(uniqueiso(beta))
	      ! Additional weighting from partitioning of elements into 'sub-types'
	      weight = weight / dsqrt((typefrac(alpha)/isofrac(uniqueiso(alpha)))*(typefrac(beta)/isofrac(uniqueiso(beta))))
	      ! Apply weighting to individual S_ab(Q)
	      do n=1,nbins
	        partialsq(alpha,beta,n) = weight * (partialsq(alpha,beta,n) - 1.0) * 4.0 * pi * numdensity
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
		  if (altbins) then
		    write(9,"(i5,E14.5)") n, partialsq(alpha,beta,n)
		  else
		    write(9,"(F6.3,3x,E14.5)") (n*binwidth)-binwidth*0.5, partialsq(alpha,beta,n)
		  end if
	        end do
	        close(9)
	      end if
	    end do
	  end do
	  write(0,*) "Finished sum."

	  ! Calculate form constant (NOT USED)
	  factor = 0.0
	  do n=1,NISOTOPES
	    if (isopops(n).eq.0) cycle
	    factor = factor + isofrac(n) * isoscatter(n)
	  end do
	  factor = factor * factor

	  ! Write total S(Q)
	  resfile=basename(1:baselen)//"sq"
	  OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	  write(9,*) "# Un-normalised data follows"
	  if (altbins) write(9,*) "# Alternative bins in use (box**2)"
	  write(9,*) "# (for information only, form constant (FZ) <b>**2 = ",factor,")"

	  do n=1,nbins
	    if (altbins) then
	      write(9,"(i5,3x,E15.7)") n,sq(n)
	    else
	      write(9,"(F15.7,E15.7)") (n-0.5)*binwidth,sq(n)
	    end if
	  end do
	  close(9)

	  ! Write elemental (isotopic) partials if required
	  if (writepartials) then
	    write(0,*) "Writing partials...."
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
		  if (altbins) then
		    write(9,"(i5,E14.5)") n, sq(n)
		  else
		    write(9,"(F6.3,3x,E14.5,I8)") (n*binwidth)-binwidth*0.5, sq(n)
		  end if
		end do
		close(9)
	      end do
	    end do
	  end if

	  close(12)
	  write(0,*) "Master finished!"
	else
	  write(0,*) "Slave ",id_mpi," finished."
	end if

999	call MPI_Finalize(err_mpi)
	end program sofq
