!	** rdfaa **
!	Calculate RDFs between specific atoms of a single species

	program rdfaa
	use parse; use dlprw; use utility
	implicit none
	integer, parameter :: maxpairs = 50
	real*8, parameter :: pi = 3.14159265358979d0
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	integer :: status,  nargs, success, baselen, npairs
	integer :: nbins,aoff1,aoff2,n,s1,s2,bin,nframes,sp,atom1(maxpairs),atom2(maxpairs)
	integer :: iargc,p,framestodo, framestodiscard = 0
	logical :: intra = .false., samemol = .false.
	real*8 :: r1(3), r2(3), r2min(3), r12(3), dist, binwidth, integral, numadded
	real*8, allocatable :: hist(:,:), rdf(:,:), sumhist(:,:), norm(:)

	binwidth=0.01   ! In Angstroms
	framestodo=-1
	sp = 1
	npairs = 0
	atom1 = 0
	atom2 = 0

	nargs = iargc()
	if (nargs.lt.5) stop "Usage : rdfaa <HISTORYfile> <OUTPUTfile> <sp> [-intra] -pair a1 a2 [-pair a1 a2 [...] ] [-frames n] [-discard n]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	sp = getargi(3)
	n = 3
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
	    case ("-pair")
	      npairs = npairs + 1
	      n = n + 1; atom1(npairs) = getargi(n)
	      n = n + 1; atom2(npairs) = getargi(n)
	    case default
	      write(0,*) "Unrecognised CLI option:", temp
	      stop
	  end select
	end do
	     
	if (npairs.eq.0) stop "No atom pairs specified!"
	do n=1,npairs
	  write(0,"(a,i2,a,2i4)") "Pair ", n, " : ", atom1(n), atom2(n)
	end do

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	nbins = maxval(cell) / binwidth
	write(0,"(A,I5)") "Target species is ",sp
	write(0,"(A,F6.3,A)") "Using binwidth of ",binwidth," Angstroms"
	write(0,"(A,I5,A)") "There will be ",nbins," histogram bins."
	if (framestodo.gt.0) write(0,"(a,i6)") "Number of frames to use in average = ",framestodo
	if (intra) write(0,"(a)") "Only intramolecular pairs will be considered"
	write(0,"(a,20(i3,'/',i3,','))") "Pairs to calculate are : ",(atom1(n),atom2(n),n=1,npairs)
	
        allocate(rdf(npairs,nbins),stat=status); if (status.GT.0) stop "Allocation error for rdf()"
	allocate(hist(npairs,nbins),stat=status); if (status.GT.0) stop "Allocation error for nrdf()"
	allocate(sumhist(npairs,nbins),stat=status); if (status.GT.0) stop "Allocation error for nrdf()"
	allocate(norm(nbins),stat=status); if (status.GT.0) stop "Allocation error for norm()"

	! Initialise the arrays...
	rdf = 0.0
	sumhist = 0

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
	numadded = 0.0
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
	
	hist = 0.0

	aoff1 = s_start(sp)
	do s1 = 1,s_nmols(sp)     ! Loop over all molecules in species

	  ! Loop over all pairs of atoms specified
	  do p = 1,npairs

	    ! Grab coordinates of what will be the central atom
	    r1(1) = xpos(aoff1+atom1(p)-1)
	    r1(2) = ypos(aoff1+atom1(p)-1)
	    r1(3) = zpos(aoff1+atom1(p)-1)
  
	    aoff2 = s_start(sp) 
	    do s2 = 1,s_nmols(sp)     ! Loop over all molecules of species 1...
	      if (s1.eq.s2) then
	        samemol = .true.
	      else
	        samemol = .false.
	      end if
  
	      ! Check for early exits
	      if ((intra.and.(.not.samemol)) .or. ((.not.intra).and.samemol)) then
	        aoff2 = aoff2 + s_natoms(sp)
	        cycle
	      end if

	      r2(1) = xpos(aoff2+atom2(p)-1)
	      r2(2) = ypos(aoff2+atom2(p)-1)
	      r2(3) = zpos(aoff2+atom2(p)-1)

	      call pbc(r2(1),r2(2),r2(3),r1(1),r1(2),r1(3),r2min(1),r2min(2),r2min(3))
	      r12 = r2min - r1
	      dist = sqrt( r12(1)*r12(1) + r12(2)*r12(2) + r12(3)*r12(3) )
  
	      bin = int(dist * (1.0/binwidth))+1
	      hist(p,bin) = hist(p,bin)+1
	      numadded = numadded+1.0
  
	      aoff2 = aoff2 + s_natoms(sp)
	    end do

	  end do ! loop over pairs

	  aoff1 = aoff1 + s_natoms(sp)

	end do   ! End main loop over all atoms of species1.

	if (nframes.EQ.1) then
	  write(0,"(A,I2,A,I4,A)") "PRDF of atoms about ",sp," : averaged over ",s_nmols(sp)," molecules."
	end if

	! Sum histogram into final RDF array
	sumhist = sumhist + hist

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
	write(0,"(a,e12.5,a,e12.5)") "Total numadded = ",numadded,", per frame = ",numadded/nframes
	write(0,"(a,e12.5)") "Total numadded per pair ",numadded / npairs

	! Normalise data
	do n=1,nbins
	  norm(n) = (4.0 * pi / 3.0) * ((n*binwidth)**3 - ((n-1)*binwidth)**3) * (s_nmols(sp) / volume(cell))
	  !norm(n) = (4.0 * pi / 3.0) * ((n*binwidth)**3 - ((n-1)*binwidth)**3) * (1.0 / volume(cell))
	  do p=1,npairs
	    rdf(p,n) = sumhist(p,n) / norm(n) / nframes / s_nmols(sp)
	  end do
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

	do p=1,npairs
	  integral = 0.0
	  resfile=basename(1:baselen)//"aardf"//CHAR(48+sp)//"_"//CHAR(48+(atom1(p)/10))//CHAR(48+MOD(atom1(p),10))// &
		& "_"//CHAR(48+(atom2(p)/10))//CHAR(48+MOD(atom2(p),10))
	  OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	  do n=1,nbins
	    integral = integral + sumhist(p,n) / nframes / s_nmols(sp)
	    write(9,"(f10.4,3x,f12.8,4x,e12.6)") (n-0.5)*binwidth,rdf(p,n),integral
	  end do
	  close(9)
	end do

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end program rdfaa

