!	** orientrdf **
!	RDF program to calculate radial distribution functions between specific orientations of molecules

	program orientrdf
	use dlprw; use utility; use parse; use IList
	implicit none
	! File names
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	! Calculation Control
	logical :: altheader = .FALSE., nonorm = .FALSE.
	integer :: framestodo = -1,frameskip = 0,framesdone, compairs(10,2), otherpairs(10,2)
	! Species Definitions
	integer :: centresp
	type(IntegerList) :: othersp
	! RDF data
	integer :: nbins
	real*8 :: binwidth,boxvolume
	real*8, allocatable :: nrdf(:,:), sumhist(:,:)
	integer, allocatable :: irdf(:,:)
	! Orientation checking
	logical :: orientcheck(MAXSP,3)			! Whether orientation check is enabled
	real*8 :: orientangle(20,3), orientdelta(20,3)	! Required orientations of molecules (if specified)
	! Working variables
	real*8, parameter :: pi = 3.14159265358979d0, radcon = 57.29577951d0
	integer :: iargc
	integer :: numadded(MAXSP),n,m,i,sp,a1,sp2,m1,m2,baselen,bin,nframes,success,o,nargs, status, sp2index
	real*8 :: dist,c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z,tx,ty,tz,bdens,integral, cellvol, norm
	real*8 :: angdelta, ax, totnumadded(MAXSP)
	logical :: failed

	binwidth=0.1   ! In Angstroms
	compairs = 0
	otherpairs = 0
	orientcheck = .FALSE.
	nargs = iargc()
	if (nargs.LT.4) stop "Usage : orientrdf <HISTORYfile> <OUTPUTfile> <sp> <othersp>[-bin width] [-header hisfile] [-frames n] [-nonorm] [-discard n] [-compair sp i j] [-otherpair sp i j] [-axis sp x1 x2 y1 y2] [-orient sp axis angle delta]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read (temp,"(i4)") centresp
	call getarg(4,temp); if (.not.parseIntegerList(temp, othersp)) stop "Failed to parse othersp list."
	
	! Write some info, probe output file, and allocate necessary arrays
	write(0,"(A,A)") "History file : ",hisfile
	write(0,"(A,A)") " Output file : ",dlpoutfile
	if (outinfo(dlpoutfile,1).eq.-1) stop "Error reading OUTPUT file."
	write(0,"(A,i3)") "Central species is = ",centresp
	write(0,"(A,20i3)") "Other species are = ", othersp%items(1:othersp%n)
	call alloc_axis()

	n = 4
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
	    case ("-axis") 
	      n = n + 1; sp = getargi(n)
	      n = n + 1; axesAatoms(sp,1) = getargi(n)
	      n = n + 1; axesAatoms(sp,2) = getargi(n)
	      n = n + 1; axesAatoms(sp,3) = getargi(n)
	      n = n + 1; axesAatoms(sp,4) = getargi(n)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp," calculated from: X=",axesAatoms(sp,1),"->", &
	        & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
	      do i=1,3
		if ((axesAatoms(sp,i).lt.1).or.(axesAatoms(sp,i).gt.s_natoms(sp))) stop "Atom id out of range for axes on this species!"
	      end do
	      axesAdefined(sp) = .true.
            case ("-bin")
              n = n + 1; binwidth = getargr(n)
              write(0,"(A,f6.2)") "Binwidth set to ",binwidth
            case ("-compair")
              n = n + 1; sp = getargi(n)
              n = n + 1; compairs(sp,1) = getargi(n)
              n = n + 1; compairs(sp,2) = getargi(n)
              write(0,"(A,3I4)") "Using COMpair for species ",sp, compairs(sp,:)
            case ("-discard")
              n = n + 1; frameskip = getargi(n)
              write(0,"(A,I4)") "Frames to discard at start: ",frameskip
            case ("-frames")
              n = n + 1; framestodo = getargi(n)
              write(0,"(A,I4)") "Frames to process: ",framestodo
            case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A,I4)") "Alternative header file supplied."
	      altheader = .TRUE.
            case ("-nonorm")
              write(0,"(A)") "RDFs will not be normalised by species number density."
	      nonorm = .TRUE.
	    case ("-orient")
	      n = n + 1; sp = getargi(n)
	      n = n + 1; m = getargi(n)
	      n = n + 1; orientangle(sp,m) = getargr(n)
	      n = n + 1; orientdelta(sp,m) = getargr(n)
	      write(0,"(A,i2,a,i1,a,f10.4,a,f10.4,a)") "Only molecules of sp ",sp," with axis ", m, " angle delta of ", &
		& orientdelta(sp,m), " degrees about ", orientangle(sp,m), " will be binned"
	      orientcheck(sp,m) = .TRUE.
            case ("-otherpair")
              n = n + 1; sp = getargi(n)
              n = n + 1; otherpairs(sp,1) = getargi(n)
              n = n + 1; otherpairs(sp,2) = getargi(n)
              write(0,"(A,3I4)") "Using other pair (for surrounding molecules) for species ",sp, otherpairs(sp,:)
	    case default
	      write(0,"(a,a)") "Unrecognised command line option:",temp
	      stop
	  end select
	end do

	! Make some sanity checks
	do sp=1,nspecies
	  ! Check that the necessary molecules have had their axes defined
	  if (axesAdefined(sp)) then
	    !write(6,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",i," calculated from: X=",axesAatoms(sp,1),"->", &
	    !  & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
	  else if (sp.eq.centresp) then
	    stop "A proper set of axes must be defined for the central species."
	  else if (orientcheck(sp,1).or.orientcheck(sp,2).or.orientcheck(sp,3)) then
	    stop "Axes must be defined on any species whose orientation is to be considered."
	  endif
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
	write(0,"(A,I5,A,F6.3,A)") "There will be ",nbins," histogram bins of ",binwidth," Angstroms."
	
	! Allocate arrays
        allocate(irdf(nspecies,nbins),stat=status); if (status.GT.0) stop "Allocation error for rdf()"
	allocate(nrdf(nspecies,nbins),stat=status); if (status.GT.0) stop "Allocation error for nrdf()"
	allocate(sumhist(nspecies,nbins),stat=status); if (status.GT.0) stop "Allocation error for sumhist()"

	! Initialise the arrays...
	nrdf = 0.0
	sumhist = 0.0
	irdf = 0

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

	! Calculate centres of mass for all molecules
	call calc_com
	do i=1,nspecies
	  ! If compairs were specified, use that instead of COM
	  if (compairs(i,1).ne.0) then
	    do m1=1,s_nmols(centresp)
	      c1x = xpos(s_start(i)+(m1-1)*s_natoms(i)+compairs(i,1)-1)
	      c1y = ypos(s_start(i)+(m1-1)*s_natoms(i)+compairs(i,1)-1)
	      c1z = zpos(s_start(i)+(m1-1)*s_natoms(i)+compairs(i,1)-1)
	      c2x = xpos(s_start(i)+(m1-1)*s_natoms(i)+compairs(i,2)-1)
	      c2y = ypos(s_start(i)+(m1-1)*s_natoms(i)+compairs(i,2)-1)
	      c2z = zpos(s_start(i)+(m1-1)*s_natoms(i)+compairs(i,2)-1)
	      call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	      comx(i,m1) = (c1x+tx)*0.5
	      comy(i,m1) = (c1y+ty)*0.5
	      comz(i,m1) = (c1z+tz)*0.5
	    end do
	  end if
	end do
	
	! Generate all molecular axes
	call genaxes()

	irdf = 0
	numadded=0

	do m1=1,s_nmols(centresp)     ! Loop over all molecules of central species

	  ! Grab centre-of-mass (or compair) coordinate for central species
	  c1x=comx(centresp,m1)
	  c1y=comy(centresp,m1)
	  c1z=comz(centresp,m1)
	
	  do sp2index=1,othersp%n

	    sp2 = othersp%items(sp2index)

	    ! Now loop over all molecules of second species....
	    do m2=1,s_nmols(sp2)

	      ! Don't add if same species *and* same molecule
	      if ((centresp.eq.sp2).and.(m1.eq.m2)) cycle

	      ! Check orientation if it has been specified
	      failed = .false.
	      do n=1,3
	        if (orientcheck(sp2,n)) then
	          ! Take dot product of axis 'n' on sp2 with that of the central molecule, and calculate the angle delta
		  if (axesBdefined(sp2)) then
		    ax = axesB(sp2,m2,(n-1)*3+1)*axesA(centresp,m1,(n-1)*3+1) + axesB(sp2,m2,(n-1)*3+2)*axesA(centresp,m1,(n-1)*3+2) + axesB(sp2,m2,(n-1)*3+3)*axesA(centresp,m1,(n-1)*3+3)
		  else
		    ax = axesA(sp2,m2,(n-1)*3+1)*axesA(centresp,m1,(n-1)*3+1) + axesA(sp2,m2,(n-1)*3+2)*axesA(centresp,m1,(n-1)*3+2) + axesA(sp2,m2,(n-1)*3+3)*axesA(centresp,m1,(n-1)*3+3)
		  end if
	          angdelta = safeAngle(ax)
	          ! If delta is negative, map delta onto +/-90
	          if ((orientdelta(sp2,n).lt.0.0).and.(angdelta.gt.90.0)) angdelta = angdelta - 180.0
	          if (dabs(angdelta-orientangle(sp2,n)).gt.dabs(orientdelta(sp2,n))) failed = .true.
		end if
		if (failed) exit
	      end do
	      if (failed) cycle

	      ! Grab the other centre coordinates if necessary
	      if (otherpairs(sp2,1).ne.0) then
	        c3x = xpos(s_start(sp2)+(m2-1)*s_natoms(sp2)+otherpairs(sp2,1)-1)
	        c3y = ypos(s_start(sp2)+(m2-1)*s_natoms(sp2)+otherpairs(sp2,1)-1)
	        c3z = zpos(s_start(sp2)+(m2-1)*s_natoms(sp2)+otherpairs(sp2,1)-1)
	        c2x = xpos(s_start(sp2)+(m2-1)*s_natoms(sp2)+otherpairs(sp2,2)-1)
	        c2y = ypos(s_start(sp2)+(m2-1)*s_natoms(sp2)+otherpairs(sp2,2)-1)
	        c2z = zpos(s_start(sp2)+(m2-1)*s_natoms(sp2)+otherpairs(sp2,2)-1)
	        call pbc(c2x,c2y,c2z,c3x,c3y,c3z,tx,ty,tz)
	        c2x = (c3x+tx)*0.5
	        c2y = (c3y+ty)*0.5
	        c2z = (c3z+tz)*0.5
	     else
		c2x=comx(sp2,m2)
		c2y=comy(sp2,m2)
		c2z=comz(sp2,m2)
	      end if

	      ! Get the shortest (MIM) distance between the two points
	      call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	      dist=sqrt( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )

	      ! 'Add' this distance...
	      bin=INT(dist/binwidth)+1
	      if (bin.lt.nbins) irdf(sp2,bin) = irdf(sp2,bin) + 1
	      numadded(sp2) = numadded(sp2) + 1
	    end do
	  end do
	end do

	! Update numadded totals
	totnumadded = totnumadded + numadded

	! Calc and store the normalised rdf....
	cellvol = volume(cell)
	do sp2=1,nspecies
	  if (nframes.EQ.1) then
	    write(0,"(A,I2,A,I2,A,I4,A)") "RDF of ",sp2," about ",centresp," : averaged over ",s_nmols(centresp)," molecules."
	  end if
	  do n=1,nbins
	    sumhist(sp2,n)=sumhist(sp2,n)+real(irdf(sp2,n))
	    if (nonorm) then
	      nrdf(sp2,n)=nrdf(sp2,n)+real(irdf(sp2,n)) / s_nmols(centresp)
	    else
	      norm = (4.0 * pi / 3.0) * ((n*binwidth)**3 - ((n-1)*binwidth)**3) * (s_nmols(sp2) / cellvol)
	      nrdf(sp2,n)=nrdf(sp2,n) + real(irdf(sp2,n)) / s_nmols(centresp) / norm
	    end if
	  end do
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

	do sp2=1,nspecies
	  if (nonorm) then
	    resfile=basename(1:baselen)//"uordf"//CHAR(48+centresp)//CHAR(48+sp2)
	  else
	    resfile=basename(1:baselen)//"ordf"//CHAR(48+centresp)//CHAR(48+sp2)
	  end if
	  open(unit=9,file=resfile,form="formatted")
	  if (compairs(centresp,1).eq.0) then
	    write(9,"(a,i2,a)") "# Central species ",centresp," calculated using all atoms for COM."
	  else
	    write(9,"(a,i2,a,i3,i3)") "# Central species ",centresp," calculated using two atoms for COM: ", compairs(centresp,1), compairs(centresp,2)
	  end if
	  if (compairs(sp2,1).eq.0) then
	    write(9,"(a,i2,a)") "# Species ",sp2," calculated using all atoms for COM."
	  else
	    write(9,"(a,i2,a,i3,i3)") "# Species ",sp2," calculated using two atoms for COM: ", compairs(sp2,1), compairs(sp2,2)
	  end if

	  integral = 0.0
	  ! Normalise the RDFs with respect to the number of frames.
	  do n=1,nbins
	    integral = integral + sumhist(sp2,n) / framesdone / s_nmols(centresp)
	    nrdf(sp2,n) = nrdf(sp2,n) / framesdone
	    write(9,"(F10.4,3(3x,F12.8))") (n*binwidth)-binwidth/2.0, nrdf(sp2,n), integral, nrdf(sp2,n)*(framesdone*s_nmols(centresp)*s_nmols(sp2)/totnumadded(sp2))
	  end do
	  close(9)
	end do
	write(0,*) "Finished."
999	close(10)
	close(13)
	end program orientrdf
