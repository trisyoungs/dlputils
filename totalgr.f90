!	** totalgr **
!	Compute the static, neutron-weighted total radial distribution function
!	Changed 13 February 2009 - Partitioning of elements into 'types' implemented

	program totgr
	use dlprw; use utility; use parse
	implicit none

	! Isotope definitions
	integer, parameter :: NISOTOPES = 11
	character*6 :: isonames(NISOTOPES)
	real*8 :: isoscatter(NISOTOPES), isofrac(NISOTOPES), isopops(NISOTOPES)
	integer :: isontypes(NISOTOPES), isonamelens(NISOTOPES)

	! DL_POLY (FF) to 'new type' to isotope mappings
	integer :: ntypes
	character*6, allocatable :: origtype(:), newtype(:), uniquetypes(:)
	integer, allocatable :: isotypes(:), typemap(:), namelens(:), uniqueiso(:)
	real*8, allocatable :: typepops(:), typefrac(:)

	! General variables
	real*8, parameter :: pi = 3.14159265358979, ftpi = 4.188790205d0
	character*80 :: hisfile,dlpoutfile,basename,resfile,namemap,headerfile
	character*20 :: temp
	integer :: i,j,k,baselen,bin,nframes,nargs,framestodo,nfftypes,alpha,beta
	integer :: n, m, o, nbins, found, frameskip, sumfac, framestodiscard, framesdone
	integer :: iargc, isuccess
	real*8 :: binwidth,boxvolume,factor,weight,delta(3),dist,svol
	logical :: success, writepartials = .FALSE., readmap = .FALSE., altheader = .FALSE.
	real*8 :: rcut,magx,magy,magz,mag,numdensity, ipos(3), jpos(3)
	real*8, allocatable :: partialgr(:,:,:),totalgr(:),weightedgr(:,:,:)

	! Scattering lengths for XX,H,D,C,N,O,F,P,S,Cl,Si
	isoscatter = (/ 0.0, -3.706, 6.671, 6.646, 9.36, 5.803, 5.654, 5.13, 2.847, 9.577, 4.1491 /)
	isonames = (/ "X ", "H ", "D ", "C ", "N ", "O ", "F ", "P ", "S ", "Cl", "Si" /)
	isonamelens = (/ 1,1,1,1,1,1,1,1,1,2,2 /)

	binwidth=0.1   ! In Angstroms
	frameskip = 1	! Take consecutive frames by default
	framestodiscard = 0

	nargs = iargc()
	if (nargs.LT.2) then
	  write(0,*) "Usage : totalgr <HISTORYfile> <OUTPUTfile> ...options"
	  write(0,*) "                                [-bin binwidth] [-frames nframes] [-discard nframes]"
	  write(0,*) "                                [-skip interval] [-partials] [-readmap <file>] [-cut cutoff]"
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
	      case ("-cut"); n = n + 1; call getarg(n,temp); read(temp,"(F20.10)") rcut
	      case ("-frames"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
	      case ("-discard"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodiscard
	      case ("-skip"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") frameskip
	      case ("-partials"); writepartials = .TRUE.
	      case ("-readmap"); readmap = .TRUE.; n = n + 1; call getarg(n,namemap)
	      case ("-altheader"); altheader = .TRUE.; n = n + 1; call getarg(n,headerfile)
	      case default
		write(0,*) "Unrecognised command line argument:", temp
	    end select
	    n = n + 1
	    if (n.GT.nargs) exit
	  end do
	end if
	
	! Open trajectory file
	call openhis(hisfile,10)
	if (readheader().EQ.-1) then
	  if (altheader) then
	    write(6,*) "Restarted trajectory:"
	    close(dlpun_his)
	    call openhis(headerfile,10)
	    if (readheader().EQ.-1) then
	      write(12,*) "Couldn't read header of alterhative history file."
	      stop
	    end if
	    close(dlpun_his)
	    call openhis(hisfile,10)
	  else
	    write(12,*) "Couldn't read header of history file."
	    goto 999
	  end if
	end if

	! CHeck imcon
	if (imcon.gt.3) stop "This image convention is not supported."

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).EQ.".") THEN
	    baselen=n
	    goto 2
	  endif
	end do
2       if (baselen.EQ.-1) THEN
	  basename="totalgr."
	  baselen=6
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=12,file=basename(1:baselen)//"grout",form="formatted",status="replace")
	write(12,"(A,F6.3,A)") "Using binwidth of ",binwidth," Angstroms"
	if (outinfo(dlpoutfile,1).EQ.-1) then
	  stop "Problem with OUTPUT file."
	end if
	write(0,*) cell

	! Read in forcefield names, the related internal types, and the corresponding isotope/element symbols
	open(unit=15,file="lengths.dat",form="formatted",status="old")
	read(15,"(I4)") nfftypes
	allocate(origtype(nfftypes))
	allocate(newtype(nfftypes))
	allocate(isotypes(nfftypes))
	allocate(typemap(natms))
	isotypes = 0
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
	    write(12,*) "Unrecognised element/isotope in lengths.dat: ", arg(3)
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
50	  write(12,"(A,I5,A,I5)") "Error reading from namefile for atom ",n," of ",natms
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

	write(0,*) cell
	rcut = cell(1) / 2.0
	nbins = rcut / binwidth + 1
	write(12,"(A,I5,A)") "There will be ",nbins," bins."

	allocate(totalgr(0:nbins))
	allocate(partialgr(ntypes,ntypes,0:nbins))
	allocate(weightedgr(ntypes,ntypes,0:nbins))

	totalgr = 0.0
	partialgr = 0.0
	weightedgr = 0.0

	! MAIN LOOP BEGINS

	nframes=0
	framesdone = 0
	numdensity = 0.0

100	isuccess=readframe()
	if (isuccess.EQ.1) then   ! End of file encountered....
	  write(12,*) "End of history file found."
	  goto 120
	end if
	if (isuccess.EQ.-1) then  ! File error....
	  write(12,*) "History file ended prematurely..."
	  goto 120
	end if
	nframes=nframes+1
	numdensity = numdensity + natms / (cell(1)*cell(5)*cell(9))
	write(0,*) nframes
	if (nframes.LE.framestodiscard) goto 100	! Discard frames at beginning of trajectory if requested
	if (mod(nframes,frameskip).ne.0) goto 100	! Frame skip interval

	! Increment partial G(r)s
	do i=1,natms-1
	  ipos(1) = xpos(i)
	  ipos(2) = ypos(i)
	  ipos(3) = zpos(i)
	  alpha = typemap(i)
	  do j=i+1,natms
	    ! Calculate MIM position of j w.r.t. i
	    call pbc(xpos(j),ypos(j),zpos(j),ipos(1),ipos(2),ipos(3),jpos(1),jpos(2),jpos(3))
	    ! Get distance
	    delta = ipos-jpos
	    dist = sqrt(delta(1)*delta(1) + delta(2)*delta(2) + delta(3)*delta(3))
	    if (dist.gt.rcut) cycle
	    bin = dist / binwidth
	    beta = typemap(j)
	    partialgr(alpha,beta,bin) = partialgr(alpha,beta,bin) + 1.0
	    partialgr(beta,alpha,bin) = partialgr(beta,alpha,bin) + 1.0
	  end do
	end do

	framesdone = framesdone + 1
	write(12,*) "Frames completed : ",framesdone

	! Next frame
	if (framesdone.eq.framestodo) goto 120
	goto 100

120	write(12,*) "Finished."
	write(12,*) ""

	! Get average number density (average over all frames read, not only those used in calc)
	numdensity = numdensity / nframes
	write(0,*) "Average total atomic number density = ",numdensity

	! Normalise partials from number density into proper RDFs
	do alpha=1,ntypes
	  do beta=1,ntypes
	    ! Number density of central atom and frame count
	write(0,"(a,2i2,2f10.4,i10)") "a/b/popa/popb/framesdone", alpha,beta,typepops(alpha), typepops(beta), framesdone
	    do n=1,nbins
	      partialgr(alpha,beta,n) =  partialgr(alpha,beta,n) / (typepops(alpha) * framesdone)
	    end do
	    ! Normalisation w.r.t. number density of second species
	    factor = typepops(beta) / (cell(1)*cell(5)*cell(9))
	    do n=1,nbins
	      svol = ftpi * ((n*binwidth)**3 - ((n-1)*binwidth)**3)
	      partialgr(alpha,beta,n) = partialgr(alpha,beta,n) / (factor*svol)
	    end do
	  end do
	end do
	write(0,*) "Finished normalisation."
	  
	! Sum into total g(r)
	totalgr = 0.0d0
	do alpha=1,ntypes
	  do beta=1,ntypes

	    ! Weight factor attributable to scattering lengths and fractional populations of species
	    weight =  typefrac(alpha)*isoscatter(uniqueiso(alpha))*typefrac(beta)*isoscatter(uniqueiso(beta))

	    ! Apply weighting to individual partials
	    do n=1,nbins
	      weightedgr(alpha,beta,n) = weight * partialgr(alpha,beta,n)
	    end do

	    ! Perform summation into totalgr
	    do n=1,nbins
	      totalgr(n) = totalgr(n) + weightedgr(alpha,beta,n)
	    end do

	    ! Write partials if required
	    if (writepartials) then
	      resfile=basename(1:baselen)//"gr"//uniquetypes(alpha)(1:namelens(alpha))//"-"//uniquetypes(beta)(1:namelens(beta))
	      OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	      do n=1,nbins
	        write(9,"(F6.3,3x,e14.5,2x,e14.5)") (n*binwidth)-binwidth*0.5, weightedgr(alpha,beta,n), partialgr(alpha,beta,n)
	      end do
	      close(9)
	    end if
	  end do
	end do
	write(0,*) "Finished summation."

	! Calculate weighting factor for G(r)
	factor = 0.0
	do alpha=1,ntypes
	  factor = factor + typefrac(alpha)*isoscatter(uniqueiso(alpha))
	end do
	factor = factor * factor

	! Write total S(Q)
	resfile=basename(1:baselen)//"totalgr"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	write(9,*) "# <b>**2 is ",factor

	do n=1,nbins
	  write(9,"(F15.7,E15.7)") (n-0.5)*binwidth,totalgr(n)/factor
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
	      totalgr = 0.0
	      do alpha=1,ntypes
	        do beta=1,ntypes
		  if (((uniqueiso(alpha).eq.i).and.(uniqueiso(beta).eq.j)).or.((uniqueiso(alpha).eq.j).and.(uniqueiso(beta).eq.i))) then
		    do n=1,nbins
		      totalgr(n) = totalgr(n) + weightedgr(alpha,beta,n)
		    end do
		  end if
		end do
	      end do
	      ! Write data
	      resfile=basename(1:baselen)//"elemgr"//isonames(i)(1:isonamelens(i))//"-"//isonames(j)(1:isonamelens(j))
	      OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	      do n=1,nbins
		write(9,"(F6.3,3x,E14.5)") (n*binwidth)-binwidth*0.5, totalgr(n)
	      end do
	      close(9)
	    end do
	  end do
	end if

999	close(12)

	end program totgr
