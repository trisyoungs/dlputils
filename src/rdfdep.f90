!	** rdfdep **
!	Calculate RDFs between specific sites of two species, dependent on a second distance

	program rdfdep
	use parse; use dlprw; use utility; use IList
	implicit none
	integer, parameter :: MAXDEPS = 20
	real*8, parameter :: pi = 3.14159265358979d0
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	integer :: status,  nargs, success, baselen
	integer :: nbins,aoff1,aoff2,n,m,m1,m2,bin,nframes,sp1,sp2,depSp
	type(IntegerList) :: sp1Site, sp2Site
	integer :: iargc, framestodo, framestodiscard = 0, nsp1used
	logical :: intra = .false., samemol = .false., logicOr = .false.
	real*8 :: i(3), j(3), jtemp(3), k(3), dist, distsq, binwidth, integral
	real*8 :: mindist(MAXDEPS), maxdist(MAXDEPS), mindistsq(MAXDEPS), maxdistsq(MAXDEPS)
	integer :: ndeps = 0, sitesOK
	type(IntegerList) :: depSite(MAXDEPS), sp1DepSite(MAXDEPS)
	integer :: nRequired(MAXDEPS), nFound(MAXDEPS)
	real*8, allocatable :: histogram(:), rdf(:), norm(:)

	binwidth=0.1   ! In Angstroms
	framestodo=-1
	nsp1used = 0

	nargs = iargc()
	if (nargs.lt.10) stop "Usage : rdfdep <HISTORYfile> <OUTPUTfile> <sp1> <atoms1> <sp2> <atoms2> <depSp> [-dep <depAtoms> <sp1atoms> <mindist> <maxdist> <n>] [-or] [-intra] [-frames n] [-discard n]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	sp1 = getargi(3)
	call getarg(4,temp)
	if (.not.parseIntegerList(temp, sp1Site)) stop "Failed to parse atom list for species 1 site."
	sp2 = getargi(5)
	call getarg(6,temp)
	if (.not.parseIntegerList(temp, sp2Site)) stop "Failed to parse atom list for species 2 site."
	depSp = getargi(7)

	n = 7
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-discard") 
	      n = n + 1; framestodiscard = getargi(n)
	    case ("-frames") 
	      n = n + 1; framestodo = getargi(n)
	    case ("-intra")
	      intra = .true.
	    case ("-dep") 
	      ndeps = ndeps + 1
	      n = n + 1; call getarg(n,temp); if (.not.parseIntegerList(temp, depSite(ndeps))) stop "Failed to parse atom list for dependent site."
	      n = n + 1; call getarg(n,temp); if (.not.parseIntegerList(temp, sp1DepSite(ndeps))) stop "Failed to parse atom list for sp1 dependent site."
	      n = n + 1; mindist(ndeps) = getargr(n)
	      mindistsq = mindist(ndeps)*mindist(ndeps)
	      n = n + 1; maxdist(ndeps) = getargr(n)
	      maxdistsq = maxdist(ndeps)*maxdist(ndeps)
	      n = n + 1; nRequired(ndeps) = getargi(n)
	    case ("-or")
	      logicOr = .true.
	    case default
	      write(0,*) "Unrecognised CLI option:", temp
	      stop
	  end select
	end do
	     
	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	nbins = maxval(cell) / binwidth
	write(0,"(A,I5,2x,a)") "Species 1 is ",sp1,s_name(sp1)
	write(0,*) "Species 1 site is (average of) : ", (sp1Site%items(n),n=1,sp1Site%n)
	write(0,"(A,I5,2x,a)") "Species 2 is ",sp2,s_name(sp2)
	write(0,*) "Species 2 site is (average of) : ", (sp2Site%items(n),n=1,sp2Site%n)
	write(0,"(A,I5,2x,a)") "Dependent species is ",depSp,s_name(depSp)
	write(0,*) "Required dependent sites are:"
	do m=1,ndeps
	  write(0,*) " ",m,". Dependent site is (average of) : ", (depSite(m)%items(n),n=1,depSite(m)%n)
	  write(0,*) "       Dependent sp1 site is (average of) : ", (sp1DepSite(m)%items(n),n=1,sp1DepSite(m)%n)
	  write(0,*) "     Dependent site min/max distance : ", mindist(m), maxdist(m)
	  write(0,*) "     Number of dependent sites required : ", nRequired(m)
	end do
	if (logicOr) write(0,*) "'Or' logic in force - any complete dependent match is enough."
	write(0,"(A,F6.3,A)") "Using binwidth of ",binwidth," Angstroms"
	write(0,"(A,I5,A)") "There will be ",nbins," histogram bins."
	if (framestodo.gt.0) write(0,"(a,i6)") "Number of frames to use in average = ",framestodo
	if (intra) write(0,"(a)") "Only intramolecular pairs will be considered"
	
        allocate(rdf(nbins),stat=status); if (status.GT.0) stop "Allocation error for rdf()"
	allocate(histogram(nbins),stat=status); if (status.GT.0) stop "Allocation error for hist()"
	allocate(norm(nbins),stat=status); if (status.GT.0) stop "Allocation error for norm()"

	! Initialise the arrays...
	rdf = 0.0
	histogram = 0

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 120  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	if (framestodiscard.gt.0) then
	  framestodiscard = framestodiscard -1
	  goto 101
	end if
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	
	aoff1 = s_start(sp1)-1
	do m1 = 1,s_nmols(sp1)     ! Loop over all molecules in species 1

	  ! Grab coordinates of the species 1 site (== i)
	  call averagePosition(sp1Site%items,sp1Site%n,aoff1,i)
  
	  ! Check for presence of dependent site(s)...
	  nFound = 0
	  sitesOK = 0
	  do m=1,ndeps
	    aoff2 = s_start(depSp)-1
	    do m2 = 1,s_nmols(depSp)

	      ! Grab coordinates of the sp1 dependent site (== k)
	      call averagePosition(sp1DepSite(m)%items,sp1DepSite(m)%n,aoff1,k)

	      ! Grab coordinates of the dependent site (== jtemp)
	      call averagePosition(depSite(m)%items,depSite(m)%n,aoff2,jtemp)
	    
	      ! Get minimum image coordinates of jtemp with k (== j)
	      call pbc(jtemp(1),jtemp(2),jtemp(3),k(1),k(2),k(3),j(1),j(2),j(3))

	      ! Get squared distance k-j
	      distsq = (k(1)-j(1))*(k(1)-j(1)) + (k(2)-j(2))*(k(2)-j(2)) + (k(3)-j(3))*(k(3)-j(3))

	      ! Is this within the specified limits
	      if ((distsq.ge.mindistsq(m)).and.(distsq.le.maxdistsq(m))) then
	        nFound(m) = nFound(m) + 1
	        if (nFound(m).ge.nRequired(m)) exit
	      end if

	      ! Increase offset counter (depSp)
	      aoff2 = aoff2 + s_natoms(depSp)

	    end do

	    ! Check here that we found all the dependent sites that were specified
	    if (nFound(m).ge.nRequired(m)) sitesOK = sitesOK + 1

	    ! Break from loop if we're done
	    if (logicOr) then
	      if (sitesOK.ge.0) exit
	    else
	      if (sitesOK.lt.m) then
	        sitesOK = 0
		exit
	      end if
	    end if

	  end do

	  ! Did we find the necessary dependent site(s) within range?
	  if (sitesOK.eq.0) then
	    aoff1 = aoff1 + s_natoms(sp1)
	    cycle
	  end if

	  ! Dependent site(s) are near enough, so accumulate RDF with species 2
	  aoff2 = s_start(sp2)-1
	  do m2 = 1,s_nmols(sp2)

	    ! Check for early exits
	    samemol = .false.
	    if (sp1.eq.sp2) samemol = (m1.eq.m2)
	    if ((intra.and.(.not.samemol)) .or. ((.not.intra).and.samemol)) then
	      aoff2 = aoff2 + s_natoms(sp2)
	      cycle
	    end if

	    ! Grab coordinates of the species 2 site (== jtemp)
	    call averagePosition(sp2Site%items,sp2Site%n,aoff2,jtemp)
	    
	    ! Get minimum image coordinates of jtemp (== j)
	    call pbc(jtemp(1),jtemp(2),jtemp(3),i(1),i(2),i(3),j(1),j(2),j(3))

	    ! Get distance i-j
	    dist = dsqrt((i(1)-j(1))*(i(1)-j(1)) + (i(2)-j(2))*(i(2)-j(2)) + (i(3)-j(3))*(i(3)-j(3)))

	    ! Add to histogram
	    bin = int(dist/binwidth)+1
	    histogram(bin) = histogram(bin)+1
  
	    ! Increase offset counter (sp2)
	    aoff2 = aoff2 + s_natoms(sp2)

	  end do

	  ! Increase sp1 counter
	  nsp1used = nsp1used + 1

	  ! Increase offset counter (sp1)
	  aoff1 = aoff1 + s_natoms(sp1)

	end do

	! Next frame
	if ((framestodo.gt.0).and.(framestodo.eq.nframes)) goto 120

	goto 101
	! Work on the results now to get the proper RDF
120	write(0,*) "Finished."
	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended prematurely!"
	goto 801
800	write(0,*) "End of unformatted HISTORY file found."

801	write(0,*) ""
	write(0,"(a,e12.5)") "Average number of sp1 selected per frame is ", real(nsp1used) / nframes

	! Normalise data
	do n=1,nbins
	  norm(n) = (4.0 * pi / 3.0) * ((n*binwidth)**3 - ((n-1)*binwidth)**3) * (s_nmols(sp2) / volume(cell))
	  !norm(n) = (4.0 * pi / 3.0) * ((n*binwidth)**3 - ((n-1)*binwidth)**3) * (1.0 / volume(cell))
	  rdf(n) = histogram(n) / norm(n) / nsp1used
	end do

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).EQ.".") THEN
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) THEN
	  basename="rdfresults."
	  baselen=11
	ELSE
	  basename=hisfile(1:baselen)
	endif

	integral = 0.0
	resfile=basename(1:baselen)//"rdfdep"//CHAR(48+sp1)//"_"//CHAR(48+sp2)
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	write(9,"(A,I5,2x,a)") "# Species 1 is ",sp1,s_name(sp1)
	write(9,*) "# Species 1 site is (average of) : ", (sp1Site%items(n),n=1,sp1Site%n)
	write(9,"(A,I5,2x,a)") "# Species 2 is ",sp2,s_name(sp2)
	write(9,*) "# Species 2 site is (average of) : ", (sp2Site%items(n),n=1,sp2Site%n)
	write(9,"(A,I5,2x,a)") "# Dependent species is ",depSp,s_name(depSp)
	write(9,*) "# Required dependent sites are:"
	do m=1,ndeps
	  write(9,*) "# ",m,". Dependent site is (average of) : ", (depSite(m)%items(n),n=1,depSite(m)%n)
	  write(9,*) "#       Dependent sp1 site is (average of) : ", (sp1DepSite(m)%items(n),n=1,sp1DepSite(m)%n)
	  write(9,*) "#     Dependent site min/max distance : ", mindist(m), maxdist(m)
	  write(9,*) "#     Number of dependent sites required : ", nRequired(m)
	end do
	if (logicOr) write(9,*) "# 'Or' logic in force - any complete dependent match is enough."
	do n=1,nbins
	  integral = integral + histogram(n) / nsp1used
	  write(9,"(f10.4,3x,f12.8,4x,e12.6,4x,e12.6)") (n-0.5)*binwidth, rdf(n), integral, histogram(n) / nsp1used
	end do
	close(9)

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end program rdfdep

