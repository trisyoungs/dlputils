!	** ardf **
!	Calculate angular RDF between sites on molecule

	program ardf
	use dlprw; use utility; use IList; use parse
	implicit none
	real*8, parameter :: pi = 3.14159265358979d0, radcon = 57.29577951d0
	character*80 :: hisfile,outfile,basename,resfile,altheaderfile
	character*20 :: temp
	logical :: altheader = .FALSE., nonorm = .FALSE., zplus = .FALSE., zminus = .FALSE.
	integer :: n,m,sp,sp1,sp2,m1,m2,baselen,bin,nframes,success,nargs, framestodo = -1,frameskip = 0,framesdone
	integer :: iargc, nbins, nangbins, angbin, aoff1, aoff2, a
	real*8 :: dist,i(3),j(3),jtemp(3), binwidth, angbinwidth, angle, ax, norm
	type(IntegerList) :: sp1Site, sp2Site
	integer, allocatable :: rdfxyz(:,:,:)

	binwidth=0.1   ! In Angstroms
	angbinwidth = 5.0   ! In degrees
	
	nargs = iargc()
	if (nargs.LT.4) stop "Usage : ardf <HISTORYfile> <OUTPUTfile> <sp1> <sp2> [-axis sp x1 x2 y1 y2] [-bin width] [-angbin width] [-header hisfile] [-frames n] [-discard n] [-sp1site atoms] [-sp2site atoms]"
	call getarg(1,hisfile)
	call getarg(2,outfile)
	sp1 = getargi(3)
	sp2 = getargi(4)
	
	if (outinfo(outfile,1).eq.-1) goto 798
	call alloc_axis()

	n = 4
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-angbin")
              n = n + 1; angbinwidth = getargr(n)
              write(0,"(A,f6.2)") "Angle binwidth set to ", angbinwidth
	    case ("-axis") 
	      n = n + 1; sp = getargi(n)
	      n = n + 1; axesAatoms(sp,1) = getargi(n)
	      n = n + 1; axesAatoms(sp,2) = getargi(n)
	      n = n + 1; axesAatoms(sp,3) = getargi(n)
	      n = n + 1; axesAatoms(sp,4) = getargi(n)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp," calculated from: X=",axesAatoms(sp,1),"->", &
	        & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
	      do m=1,3
		if ((axesAatoms(sp,m).lt.1).or.(axesAatoms(sp,m).gt.s_natoms(sp))) stop "Atom id out of range for axes on this species!"
	      end do
	      axesAdefined(sp) = .true.
            case ("-bin")
              n = n + 1; binwidth = getargr(n)
              write(0,"(A,f6.2)") "Binwidth set to ",binwidth
            case ("-discard")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") frameskip
              write(0,"(A,I4)") "Frames to discard at start: ",frameskip
            case ("-frames")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
              write(0,"(A,I4)") "Frames to process: ",framestodo
            case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A,I4)") "Alternative header file supplied."
	      altheader = .TRUE.
	    case ("-sp1site")
	      n = n + 1; sp = getargi(n)
	      if (sp1Site%n.gt.0) stop "Definition of site specified twice for species 1"
	      n = n + 1; call getarg(n,temp); if (.not.parseIntegerList(temp, sp1Site)) stop "Failed to parse atom list."
	      write(0,"(i2,a,i2,a,20i4)") sp1Site%n, " atoms added as other sites for species ", sp, ": ", sp1Site%items(1:sp1Site%n)
	    case ("-sp2site")
	      n = n + 1; sp = getargi(n)
	      if (sp2Site%n.gt.0) stop "Definition of site specified twice for species 2"
	      n = n + 1; call getarg(n,temp); if (.not.parseIntegerList(temp, sp2Site)) stop "Failed to parse atom list."
	      write(0,"(i2,a,i2,a,20i4)") sp2Site%n, " atoms added as other sites for species ", sp, ": ", sp2Site%items(1:sp2Site%n)
	    case default
	      write(0,"(a,a)") "Unrecognised command line option:",temp
	      stop
	  end select
	end do

        ! Now, read in the history header so that we have cell()
	! If this fails then we may have a restarted trajectory. Continue, but only if
	! a header can be read in from the specified alternative history file..
	call openhis(hisfile,10)
        if (readheader().EQ.-1) then
	  if (altheader) then
	    write(0,*) "Restarted trajectory:"
	    close(dlpun_his)
	    call openhis(altheaderfile,10)
	    if (readheader().EQ.-1) goto 797
	    close(dlpun_his)
	    call openhis(hisfile,10)
	  else
	    goto 797
	  end if
	end if

	nbins = maxval(cell) / binwidth + 1
	nangbins = 180.0 / angbinwidth + 1
	write(0,"(A,I5,A,F6.3,A)") "There will be ",nbins," distance histogram bins of ",binwidth," Angstroms."
	write(0,"(A,I5,A,F6.3,A)") "There will be ",nangbins," angle histogram bins of ",angbinwidth," Degrees."
	allocate(rdfxyz(3,nbins,nangbins))
	rdfxyz = 0

	! Check that the relevant axes definitions have been supplied
	do sp=1,nspecies
	  ! Check that the necessary molecules have had their axes defined
	  if (axesAdefined(sp)) then
	    write(6,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp," calculated from: X=",axesAatoms(sp,1),"->", &
	      & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
	  else if ((sp.eq.sp1).or.(sp.eq.sp2)) then
	    stop "A proper set of axes must be defined for both sp1 and sp2."
	  end if
	end do

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
100	nframes=0
	framesdone = 0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes, framesdone
	if (nframes.le.frameskip) goto 101

	framesdone = framesdone + 1

	! Generate all molecular axes
	call genaxes()

	aoff1 = s_start(sp1) - 1
	do m1=1,s_nmols(sp1)     ! Loop over all molecules of species 1...

	  ! Grab coordinate for central species
	  if (sp1Site%n.eq.0) then
	    i(1)=axesAorigin(sp1,m1,1)
	    i(2)=axesAorigin(sp1,m1,2)
	    i(3)=axesAorigin(sp1,m1,3)
	  else
	    call averagePosition(sp1Site%items,sp1Site%n,aoff1,i)
	  end if
	
	  ! Now loop over all molecules of second species....
	  aoff2 = s_start(sp2) - 1
	  do m2=1,s_nmols(sp2)

	    ! Don't add if same species *and* same molecule
	    if ((sp1.eq.sp2).and.(m1.eq.m2)) then
	      aoff2 = aoff2 + s_natoms(sp2)
	      cycle
	    end if

	    ! Grab the other coordinates
	    if (sp1Site%n.eq.0) then
	      jtemp(1)=axesAorigin(sp2,m2,1)
	      jtemp(2)=axesAorigin(sp2,m2,2)
	      jtemp(3)=axesAorigin(sp2,m2,3)
	    else
	      call averagePosition(sp2Site%items,sp2Site%n,aoff2,jtemp)
	    end if
	    
	    ! Get minimum image coordinates of j, and get distance i-j
	    call pbc(jtemp(1),jtemp(2),jtemp(3),i(1),i(2),i(3),j(1),j(2),j(3))
	    dist = sqrt( (i(1)-j(1))*(i(1)-j(1)) + (i(2)-j(2))*(i(2)-j(2)) + (i(3)-j(3))*(i(3)-j(3)))
	    bin = int(dist/binwidth)+1
	    if (bin.gt.nbins) then
	      aoff2 = aoff2 + s_natoms(sp2)
	      cycle
	    end if

	    ! Calculate axis angles
	    do n=1,3
	      ! Take dot product of axis on sp2 with that of the central molecule, and calculate the angle delta
	      ax = axesA(sp2,m2,(n-1)*3+1)*axesA(sp1,m1,(n-1)*3+1) + axesA(sp2,m2,(n-1)*3+2)*axesA(sp1,m1,(n-1)*3+2) + axesA(sp2,m2,(n-1)*3+3)*axesA(sp1,m1,(n-1)*3+3)
	      angle = acos(ax) * radcon
	      angbin = int(angle/angbinwidth)+1
	      rdfxyz(n,bin,angbin) = rdfxyz(n,bin,angbin) + 1
	    end do

	    ! Increase atom offset
	    aoff2 = aoff2 + s_natoms(sp2)

	  end do

	  ! Increase atom offset
	  aoff1 = aoff1 + s_natoms(sp1)

	end do

	if (framesdone.eq.framestodo) goto 801
	! Next frame
	goto 101

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

797	write(0,*) "No header found in history file. If a restarted trajectory, use '-header'"
	goto 999
798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended prematurely!"
	goto 801
800	write(0,*) "End of unformatted HISTORY file found."
801	write(0,*) ""

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="rdfresults."
	  baselen=11
	else
	  basename=hisfile(1:baselen)
	endif

	! Loop over axes
	do a=1,3
	  resfile=basename(1:baselen)//"ardf"//CHAR(48+sp1)//CHAR(48+sp2)//CHAR(119+a)
	  open(unit=9,file=resfile,form="formatted")
	  if (sp1Site%n.eq.0) then
	    write(9,"(a,i2,a)") "# Species ",sp1," site is axis origin."
	  else
	    write(0,"(a,i2,a,20i4)") "# Species ", sp1, ", site is average of: ", sp1Site%items(1:sp1Site%n)
	  end if
	  if (sp2Site%n.eq.0) then
	    write(9,"(a,i2,a)") "# Species ",sp2," site is axis origin."
	  else
	    write(0,"(a,i2,a,20i4)") "# Species ", sp2, ", site is average of: ", sp2Site%items(1:sp2Site%n)
	  end if

	  ! Normalise the RDFs with respect to the number of frames and solid angle.
	  do m=1,nangbins
	    do n=1,nbins
	      norm = (4.0/3.0) * pi * ((n*binwidth)**3 - ((n-1)*binwidth)**3) * s_nmols(sp1) / volume(cell)
	      write(9,"(F10.4,2(3x,F12.8))") ((n-0.5)*binwidth), rdfxyz(a,n,m) / norm / nframes / s_nmols(sp1) / dsin((m-0.5)*angbinwidth/RADCON)
	    end do
	    write(9,*) ""
	  end do
	  close(9)
	end do

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program ardf
