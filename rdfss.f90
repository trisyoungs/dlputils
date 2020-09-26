!	** rdfss **
!	RDF program to calculate pair-pair distribution functions between pairs of atoms
!	or a single atom and the gemetric / mass centre of a second species.
!	Calculates the density at each frame of the simulation (for NPT simulations)
!	Default is that all site-site g(r)'s between the sp1 atoms and the COM of the
!	second species is required.

	program rdfss
	use parse; use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	integer :: n,a1,s1,m1,m2,baselen,bin,nframes,success,o,nargs,numadded,sp1,sp2,compair(2)
	integer :: framestodiscard = 0, framesdone = 0
	integer :: iargc, status
	real*8 :: dist,c1x,c1y,c1z,c2x,c2y,c2z,tx,ty,tz,numdens,const,norm
	integer :: nbins
	real*8 :: binwidth,boxvolume
	real*8, allocatable :: rdf(:,:),purerdf(:),nrdf(:,:)

	binwidth=0.1   ! In Angstroms
	sp1 = 1
	sp2 = 2
	compair = 0

	nargs = iargc()
	if (nargs.LT.2) stop "Usage : rdfss <HISTORYfile> <OUTPUTfile> <sp1> <sp2> [-bin binwidth] [-compair i j] [-discard n]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	sp1 = getargi(3)
	sp2 = getargi(4)
	if (nargs.GE.5) then
	  n = 5
	  do
	    call getarg(n,temp)
	    select case (temp)
	      case ("-bin"); n = n + 1; binwidth = getargr(n)
              case ("-compair")
                n = n + 1; compair(1) = getargi(n)
                n = n + 1; compair(2) = getargi(n)
                write(0,"(A,3I4)") "Using COMpair for species 2 ",compair(:)
              case ("-discard")
                n = n + 1; framestodiscard = getargi(n)
                write(0,"(a,i6)") "Frames to discard at start of trajectory:", framestodiscard
	      case default
		write(0,*) "Unrecognised argument: ", temp
		stop
	    end select
	    n = n + 1
	    if (n.GT.nargs) exit
	  end do
	end if
	
	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	nbins = cell(1) / binwidth
	write(0,"(A,F6.3,A)") "Using binwidth of ",binwidth," Angstroms"
	write(0,"(A,I5,A)") "There will be ",nbins," histogram bins."
	write(0,"(A,I1,A,I1)") "Calculating PRDFs between atoms of species ",sp1," and centre-of-mass of species ",sp2
	
	! Initialise the arrays...
        allocate(rdf(maxval(s_natoms),nbins),stat=status); if (status.GT.0) stop "Allocation error for rdf()"
        allocate(purerdf(nbins),stat=status); if (status.GT.0) stop "Allocation error for purerdf()"
	allocate(nrdf(maxval(s_natoms),nbins),stat=status); if (status.GT.0) stop "Allocation error for nrdf()"
	rdf = 0.0
	purerdf = 0.0
	nrdf = 0.0

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 800  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	if (nframes.le.framestodiscard) goto 101
	
	call calc_com
	if (compair(1).ne.0) then
	  do m2=1,s_nmols(sp2)
	    c1x = xpos(s_start(sp2)+(m2-1)*s_natoms(sp2)+compair(1)-1)
	    c1y = ypos(s_start(sp2)+(m2-1)*s_natoms(sp2)+compair(1)-1)
	    c1z = zpos(s_start(sp2)+(m2-1)*s_natoms(sp2)+compair(1)-1)
	    c2x = xpos(s_start(sp2)+(m2-1)*s_natoms(sp2)+compair(2)-1)
	    c2y = ypos(s_start(sp2)+(m2-1)*s_natoms(sp2)+compair(2)-1)
	    c2z = zpos(s_start(sp2)+(m2-1)*s_natoms(sp2)+compair(2)-1)
	    call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	    comx(sp2,m2) = (c1x+tx)*0.5
	    comy(sp2,m2) = (c1y+ty)*0.5
	    comz(sp2,m2) = (c1z+tz)*0.5
	  end do
	end if

	do a1=1,s_natoms(sp1)     ! Loop over all atoms of species 1....
	  rdf(a1,1:nbins) = (/ (0,n=1,nbins) /)
	  o=s_start(sp1)	 ! Set the start atom of the species
	  do m1=1,s_nmols(sp1)     ! Loop over all molecules of species 1...
	    c1x=xpos(o+a1-1)
	    c1y=ypos(o+a1-1)
	    c1z=zpos(o+a1-1)
	    numadded=0
	    ! Now loop over all molecules of second species....
	    do m2=1,s_nmols(sp2)
	      ! Set the coordinates for the second point....
	      ! Grab the GC value grom the GCxyz arrays....
	      c2x=comx(sp2,m2)
	      c2y=comy(sp2,m2)
	      c2z=comz(sp2,m2)
	      ! Get the shortest (MIM) distance between the atom pair...
	      call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	      dist=sqrt( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )
	      ! 'Add' this distance...
	      ! bin=INT(dist*10)+1
	      bin = INT(dist * (1.0/binwidth))+1
	      rdf(a1,bin)=rdf(a1,bin)+1
	      numadded=numadded+1
	    end do
            ! Consistency check - we should *always* add s_nmols() unless s1 = s2.
            if (sp1.eq.sp2) then
!              if (numadded.NE.s_nmols(sp1)-1) write(0,*) &
 !               & "WARNING: Did not bin all molecules...",nframes,sp1,sp2,numadded
            else
              if (numadded.NE.s_nmols(sp2)) write(0,*) &
                & "WARNING: Did not bin all molecules...",nframes,sp1,sp2,numadded
            end if
	    o=o+s_natoms(sp1)
	  end do
	end do   ! End main loop over all atoms of sp1.

	if (framesdone.EQ.1) write(0,*) "numadded:=",numadded
	if (framesdone.EQ.1) then
	  if (sp1.EQ.sp2) write(0,"(A,I2,A,I4,A)") "PRDF of atoms about ",sp2," : averaged over ",s_nmols(sp2)-1," molecules."
	  if (sp1.NE.sp2) write(0,"(A,I2,A,I4,A)") "PRDF of atoms about ",sp2," : averaged over ",s_nmols(sp2)," molecules."
	end if
	if (sp1.EQ.sp2) rdf = rdf / (s_nmols(sp2)-1)
	if (sp1.NE.sp2) rdf = rdf / s_nmols(sp2)
	nrdf = nrdf + rdf

	! Next frame
	framesdone = framesdone + 1
	goto 101

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
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

	numdens = s_nmols(sp1) / volume(cell)
	const = (4.0*3.141592654) / 3.0
	do a1=1,s_natoms(sp1)
	  resfile=basename(1:baselen)//"ssrdf"//CHAR(48+sp2)//"_"//CHAR(48+sp1)//"_"//CHAR(48+(a1/10))//CHAR(48+MOD(a1,10))//atmname(s_start(sp1)+a1-1)(1:2)
	  OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	  !write(9,*) "Bin     Raw     G(r)"
	  ! Normalise the RDFs with respect to the number of frames and number density of species
	  do n=1,nbins
	    norm = (const*((n*binwidth)**3-((n-1)*binwidth)**3)) * numdens
	    nrdf(a1,n)=nrdf(a1,n)/framesdone / norm
	  end do
	  do n=1,nbins
	    ! write(9,"(F6.3,3x,2F12.8)") (n*binwidth)-binwidth/2.0,purerdf(n),nrdf(a1,n)
	    write(9,"(F6.3,3x,F12.8)") (n*binwidth)-binwidth/2.0,nrdf(a1,n)
	  end do
	  close(9)
	end do
	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end program rdfss

