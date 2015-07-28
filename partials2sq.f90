!	** raw2sq **
!	Compute the static, neutron-weighted Faber-Ziman total structure factor S(Q) from raw, unweighted partials

	program raw2sq
	use dlprw; use utility; use parse
	implicit none

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
	character*80 :: dlpoutfile,basename,resfile,namemap,hisfile,headerfile,lengthsfile
	character*20 :: temp
	integer :: i,j,k,baselen,bin,nargs,nfftypes,alpha,beta
	integer :: n, m, o, nbins, found
	integer :: iargc
	logical :: success, writepartials = .FALSE., readmap = .FALSE., altheader = .FALSE.
	real*8 :: kcut, binwidth, factor, selfscatter
	real*8, allocatable :: partialsq(:,:,:),sq(:)

	! Scattering lengths for H,D,C,N,O,F,P,S,Cl
	isoscatter = (/ -3.706, 6.671, 6.646, 9.36, 5.803, 5.654, 5.13, 2.847, 9.577 /)
	isonames = (/ "H ", "D ", "C ", "N ", "O ", "F ", "P ", "S ", "Cl" /)
	isonamelens = (/ 1,1,1,1,1,1,1,1,2 /)

	binwidth=0.1   ! In Angstroms
	kcut = 5.0    ! Reciprocal space cutoff (box integers)
        lengthsfile="lengths.dat"

	nargs = iargc()
	if (nargs.LT.2) then
	  write(0,*) "Usage : raw2sq <HISTORYfile> <OUTPUTfile> ...options"
	  write(0,*) "            [-bin width]        Set x binwidth to use in S(Q) output"
	  write(0,*) "            [-kcut cutoff]      Set maximum k magnitude to bin"
	  write(0,*) "            [-partials]         Write partial S(Q) info"
	  write(0,*) "            [-readmap <file>]   Read alternative atom names map from <file>"
	  write(0,*) "            [-altheader <file>] Use specified DL_POLY history <file> for header"
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
	      case ("-readmap"); readmap = .TRUE.; n = n + 1; call getarg(n,namemap)
	      case ("-altheader"); altheader = .TRUE.; n = n + 1; call getarg(n,headerfile)
	      case ("-partials"); writepartials = .TRUE.
	      case ("-lengths"); n = n + 1; call getarg(n,lengthsfile)
	      case default
		write(0,*) "Unrecognised command line argument:", temp
		stop
	    end select
	    n = n + 1
	    if (n.GT.nargs) exit
	  end do
	end if
	
	nbins = kcut / binwidth + 1
 
	! Open trajectory file
	call openhis(hisfile,10)
	if (readheader().EQ.-1) then
	  if (altheader) then
	    write(0,*) "Restarted trajectory:"
	    close(dlpun_his)
	    call openhis(headerfile,10)
	    if (readheader().EQ.-1) then
	      write(0,*) "Couldn't read header of alterhative history file."
	      goto 999
	    end if
	    close(dlpun_his)
	    call openhis(hisfile,10)
	  else
	    write(0,*) "Couldn't read header of history file."
	    goto 999
	  end if
	end if

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).EQ.".") THEN
	    baselen=n
	    goto 2
	  endif
	end do
2       if (baselen.EQ.-1) THEN
	  basename="sq."
	  baselen=6
	else
	  basename=hisfile(1:baselen)
	endif

	write(0,"(A,F6.3,A)") "Using binwidth of ",binwidth," Angstroms"
	write(0,"(A,I5,A)") "There will be ",nbins," histogram bins."
	if (outinfo(dlpoutfile,1).EQ.-1) then
	  write(0,*) "Problem with OUTPUT file."
	  goto 999
	end if

	! Read in forcefield names, the related internal types, and the corresponding isotope/element symbols
	open(unit=15,file=lengthsfile,form="formatted",status="old")
	read(15,"(I4)") nfftypes
	allocate(origtype(nfftypes))
	allocate(newtype(nfftypes))
	allocate(isotypes(nfftypes))
	allocate(typemap(natms))
	! First column is original atom name, second is new type name, third is related element/isotope
	do n=1,nfftypes
	  success = readline(15)
	  origtype(n) = arg(1)
	  newtype(n) = arg(2)
	  !read(15,"(A6,A6,)") origtype(n),newtype(n),temp
	  isotypes(n) = 0
	  do m=1,NISOTOPES
	    if (arg(3).eq.isonames(m)) isotypes(n) = m
	  end do
	  if (isotypes(n).eq.0) then
	    write(0,*) "Unrecognised element/isotope in lengths.dat: ",temp
	    stop
	  end if
	end do
	close(15)

	! If specified, read in a different set of atom names from a specified file (-readmap)
	! Allows for the changing of atom names or easier splitting of FF types to 1 or more types
	if (readmap) then
	  write(0,"(A,A)") "Reading alternative list of atom names from: ",namemap
	  open(unit=15,file=namemap,form="formatted",status="old",err=50)
	  do n=1,natms
	    read(15,"(A6)",err=50) atmname(n)
	  end do
	  goto 55
50	  write(0,"(A,I5,A,I5)") "Error reading from namefile for atom ",n," of ",natms
55	  close(15)
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
	    write(0,*) "Unrecognised element/isotope in unique list - this shouldn't have happened! : ",uniquetypes(n)
	    stop
	  end if
	  isontypes(uniqueiso(n)) = isontypes(uniqueiso(n)) + 1
	end do
	write(0,"(A,I2,A)") "Atoms are to be partitioned into ",ntypes," unique groups:"
	  
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

	write(0,*) "   Type      N      Frac   El/Iso   S.Length"
	write(0,"(i2,2x,a6,1x,f7.0,3x,f7.4,2x,a6,2x,f7.4)") (n,uniquetypes(n),typepops(n),typepops(n)/natms,isonames(uniqueiso(n)),isoscatter(uniqueiso(n)),n=1,ntypes)

	write(0,*) ""
	write(0,*) "Element/isotope fractional populations:"
	write(0,*) "  Name     Frac    S.Length   NTypes"
	do n=1,NISOTOPES
	  if (isopops(n).eq.0) cycle
	  write(0,"(i2,2x,a6,1x,f7.5,3x,f7.4,2x,i4)") n,isonames(n),isofrac(n),isoscatter(n),isontypes(n)
	end do

	! Set the name lengths of each of the typenames
	do n=1,ntypes
	  do m=6,1,-1
	    if (uniquetypes(n)(m:m).EQ." ") namelens(n) = m-1
	  end do
	end do

	deallocate(newtype)
	deallocate(origtype)

	allocate(sq(nbins))
	allocate(partialsq(ntypes,ntypes,nbins))

	partialsq = 0.0
	sq = 0.0d0
	selfscatter = 0.0

	! Sum into total S(Q)
	do alpha=1,ntypes
	  do beta=1,ntypes

	    ! Load raw data from file and weight by scattering length and new atomic composition
	    resfile=basename(1:baselen)//"partsq"//uniquetypes(alpha)(1:namelens(alpha))//"-"//uniquetypes(beta)(1:namelens(beta))
	    open(unit=15,form="formatted",status="old",file=resfile,ERR=999)
	    ! Skip comments at start of file
	    do n=1,3
	      if (.not.readline(15)) goto 999
	    end do
	    do n=1,nbins
	      if (.not.readline(15)) goto 999
	      partialsq(alpha,beta,n) = argr(2)
	      !read(15,"(9x,f20.14)") partialsq(alpha,beta,n)
	      !partialsq(alpha,beta,n) = weight * partialsq(alpha,beta,n)
	    end do
	    close(15)

	    ! Calculate weighting factor for partial (does not include fractional population since this is accounted for)
	    factor = isoscatter(uniqueiso(alpha)) * isoscatter(uniqueiso(beta))

	    ! Increment self-scattering correction
	    if (alpha.eq.beta) selfscatter = selfscatter + typefrac(alpha) * factor

	    ! Perform summation into total S(Q)
	    write(6,"(a,f9.6,a,f9.6,a,f12.6)") "Partial "//uniquetypes(alpha)(1:namelens(alpha))//"-"//uniquetypes(beta)(1:namelens(beta))//" has fractional pops ",typefrac(alpha)," and ",typefrac(beta), "for a total atomic weighting of ", isoscatter(uniqueiso(alpha))*isoscatter(uniqueiso(beta)) / 100.0
	    do n=1,nbins
	      sq(n) = sq(n) + partialsq(alpha,beta,n) * factor
	    end do
	  end do
	end do
	write(0,*) "Finished sum"

	! Remove self-scattering from S(Q) and convert between units of fm (www.ncnr.nist.gov) and 
	! 10**-12 cm units (ISIS) used for the scattering lengths
	sq = (sq - selfscatter) / 100.0

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
	write(9,*) "# Un-normalised data follows."
	write(9,*) "# (for information only, form constant (FZ) <b>**2 = ",factor,")"
	write(9,*) "#   Type      N      Frac   El/Iso   S.Length"
	write(9,"('# ',i2,2x,a6,1x,f7.0,3x,f7.4,2x,a6,2x,f7.4)") (n,uniquetypes(n),typepops(n),typepops(n)/natms,isonames(uniqueiso(n)),isoscatter(uniqueiso(n)),n=1,ntypes)
	write(9,*) "# Element/isotope fractional populations:"
	write(9,*) "#  Name     Frac    S.Length   NTypes"
	do n=1,NISOTOPES
	  if (isopops(n).eq.0) cycle
	  write(9,"('# ',i2,2x,a6,1x,f7.5,3x,f7.4,2x,i4)") n,isonames(n),isofrac(n),isoscatter(n),isontypes(n)
	end do

	do n=1,nbins
	  write(9,"(F15.7,E15.7,I8)") (n-0.5)*binwidth,sq(n)
	end do
	close(9)

	stop
999	write(0,"(a,a)") "Error reading from file: ",resfile
	end program raw2sq
