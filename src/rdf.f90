!	** rdf **
!	RDF program to calculate radial distribution functions between the centres-of-mass of
!	molecules. Automatically performed for all species pair combinations.
!	Calculates the density at each frame of the simulation (for NPT simulations)

	program rdf
	use dlprw; use utility; use parse
	implicit none
	! File names
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	! Calculation Control
	logical :: altheader = .FALSE., nonorm = .FALSE., zplus = .FALSE., zminus = .FALSE.
	integer :: framestodo = -1,frameskip = 0,framesdone, compairs(10,2), otherpairs(10,2)
	! RDF data
	integer :: nbins
	real*8 :: binwidth,boxvolume
	real*8, allocatable :: purerdf(:),nrdf(:,:,:),sumhist(:,:,:)
	integer, allocatable :: irdf(:,:,:)
	! Working variables
	real*8, parameter :: pi = 3.14159265358979d0
	integer :: iargc
	integer :: numadded,n,m,i,a1,sp1,sp2,m1,m2,baselen,bin,nframes,success,o,nargs, status
	real*8 :: dist,c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z,tx,ty,tz,bdens,integral, cellvol, norm

	binwidth=0.1   ! In Angstroms
	compairs = 0
	otherpairs = 0
	nargs = iargc()
	if (nargs.LT.2) stop "Usage : rdf <HISTORYfile> <OUTPUTfile> [-bin width] [-header hisfile] [-frames n] [-nonorm] [-zplus] [-zminus] [-discard n] [-compair sp i j] [-otherpair sp i j]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	
	n = 2
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-bin")
              n = n + 1; binwidth = getargr(n)
              write(0,"(A,f6.2)") "Binwidth set to ",binwidth
            case ("-compair")
              n = n + 1; sp1 = getargi(n)
              n = n + 1; compairs(sp1,1) = getargi(n)
              n = n + 1; compairs(sp1,2) = getargi(n)
              write(0,"(A,3I4)") "Using COMpair for species ",sp1, compairs(sp1,:)
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
            case ("-otherpair")
              n = n + 1; sp1 = getargi(n)
              n = n + 1; otherpairs(sp1,1) = getargi(n)
              n = n + 1; otherpairs(sp1,2) = getargi(n)
              write(0,"(A,3I4)") "Using other pair (for surrounding molecules) for species ",sp1, otherpairs(sp1,:)
            case ("-zminus")
              write(0,"(A)") "Minimum image molecules will only be accepted if z < zcentre."
	      if (zplus) stop "Can't use both -zplus and -zminus"
	      zminus = .TRUE.
            case ("-zplus")
              write(0,"(A)") "Minimum image molecules will only be accepted if z > zcentre."
	      if (zminus) stop "Can't use both -zplus and -zminus"
	      zplus = .TRUE.
	    case default
	      write(0,"(a,a)") "Unrecognised command line option:",temp
	      stop
	  end select
	end do

	! Open and check the files...
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798

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
        allocate(irdf(nspecies,nspecies,nbins),stat=status); if (status.GT.0) stop "Allocation error for rdf()"
        allocate(purerdf(nbins),stat=status); if (status.GT.0) stop "Allocation error for purerdf()"
	allocate(nrdf(nspecies,nspecies,nbins),stat=status); if (status.GT.0) stop "Allocation error for nrdf()"
	allocate(sumhist(nspecies,nspecies,nbins),stat=status); if (status.GT.0) stop "Allocation error for sumhist()"

	! Initialise the arrays...
	nrdf = 0.0
	sumhist = 0.0
	irdf = 0
	purerdf = 0.0

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
	do sp1=1,nspecies
	  ! If compairs were specified, use that instead of COM
	  if (compairs(sp1,1).ne.0) then
	    do m1=1,s_nmols(sp1)
	      c1x = xpos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,1)-1)
	      c1y = ypos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,1)-1)
	      c1z = zpos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,1)-1)
	      c2x = xpos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,2)-1)
	      c2y = ypos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,2)-1)
	      c2z = zpos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,2)-1)
	      call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	      comx(sp1,m1) = (c1x+tx)*0.5
	      comy(sp1,m1) = (c1y+ty)*0.5
	      comz(sp1,m1) = (c1z+tz)*0.5
	    end do
	  end if
	end do
	
	irdf = 0
	do sp1=1,nspecies

	  do m1=1,s_nmols(sp1)     ! Loop over all molecules of species 1...

	    ! Grab centre-of-mass (or compair) coordinate for central species
	    c1x=comx(sp1,m1)
	    c1y=comy(sp1,m1)
	    c1z=comz(sp1,m1)
	
	    do sp2=1,nspecies
	      numadded=0
	      ! Now loop over all molecules of second species....
	      do m2=1,s_nmols(sp2)
		! Don't add if same species *and* same molecule
		if ((sp1.eq.sp2).and.(m1.eq.m2)) cycle
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
		! Check zplus or zminus
		if (zplus.and.(tz.gt.c1z)) cycle
		if (zminus.and.(tz.lt.c1z)) cycle
		dist=sqrt( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )
		! 'Add' this distance...
		bin=INT(dist/binwidth)+1
		if (bin.lt.nbins) irdf(sp1,sp2,bin) = irdf(sp1,sp2,bin) + 1
		numadded=numadded+1
	      end do
	      ! Consistency check - we should *always* add s_nmols() unless sp1 = sp2.
	      ! Don't do this check if zplus or zminus is in effect
	      if ((.not.zplus).and.(.not.zminus)) then
	        if (sp1.eq.sp2) then
	          if (numadded.NE.s_nmols(sp1)-1) write(0,"(A,i6,a,i2,a,i2,' (',i7,')')") &
		    & "WARNING: Frame ",nframes," - Did not bin all molecules of ",sp2," around ",sp1,numadded
	        else
		  if (numadded.NE.s_nmols(sp2)) write(0,*) &
		    & "WARNING: Did not bin all molecules...",nframes,sp2,numadded
		end if
	      end if
	    end do
	  end do
	end do

	! Calc and store the normalised rdf....
	cellvol = volume(cell)
	do sp2=1,nspecies
	  do sp1=1,nspecies
	    if (nframes.EQ.1) then
	      write(0,"(A,I2,A,I2,A,I4,A)") "RDF of ",sp2," about ",sp1," : averaged over ",s_nmols(sp1)," molecules."
	    end if
	    do n=1,nbins
	      sumhist(sp1,sp2,n)=sumhist(sp1,sp2,n)+real(irdf(sp1,sp2,n))
	      if (nonorm) then
	        nrdf(sp1,sp2,n)=nrdf(sp1,sp2,n)+real(irdf(sp1,sp2,n)) / s_nmols(sp1)
	      else
	        norm = (4.0 * pi / 3.0) * ((n*binwidth)**3 - ((n-1)*binwidth)**3) * (s_nmols(sp2) / cellvol)
	        nrdf(sp1,sp2,n)=nrdf(sp1,sp2,n)+real(irdf(sp1,sp2,n)) / s_nmols(sp1) / norm
	      end if
	    end do
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

	do sp1=1,nspecies
	  do sp2=1,nspecies
	    if (nonorm) then
	      if (zplus) then
		resfile=basename(1:baselen)//"u_zp_rdf"//CHAR(48+sp1)//CHAR(48+sp2)
	      else if (zminus) then
		resfile=basename(1:baselen)//"u_zm_rdf"//CHAR(48+sp1)//CHAR(48+sp2)
	      else
		resfile=basename(1:baselen)//"urdf"//CHAR(48+sp1)//CHAR(48+sp2)
	      end if
	    else
	      if (zplus) then
		resfile=basename(1:baselen)//"zp_rdf"//CHAR(48+sp1)//CHAR(48+sp2)
	      else if (zminus) then
		resfile=basename(1:baselen)//"zm_rdf"//CHAR(48+sp1)//CHAR(48+sp2)
	      else
		resfile=basename(1:baselen)//"rdf"//CHAR(48+sp1)//CHAR(48+sp2)
	      end if
	    end if
	    open(unit=9,file=resfile,form="formatted")
	    if (compairs(sp1,1).eq.0) then
	      write(9,"(a,i2,a)") "# Species ",sp1," calculated using all atoms for COM."
	    else
	      write(9,"(a,i2,a,i3,i3)") "# Species ",sp1," calculated using two atoms for COM: ", compairs(sp1,1), compairs(sp1,2)
	    end if
	    if (compairs(sp2,1).eq.0) then
	      write(9,"(a,i2,a)") "# Species ",sp2," calculated using all atoms for COM."
	    else
	      write(9,"(a,i2,a,i3,i3)") "# Species ",sp2," calculated using two atoms for COM: ", compairs(sp2,1), compairs(sp2,2)
	    end if

	    integral = 0.0
	    ! Normalise the RDFs with respect to the number of frames.
	    do n=1,nbins
	      integral = integral + sumhist(sp1,sp2,n) / framesdone / s_nmols(sp1)
	      nrdf(sp1,sp2,n) = nrdf(sp1,sp2,n) / framesdone
	      ! write(9,"(F6.3,3x,F12.8)") (n*binwidth)-binwidth/2.0,nrdf(a1,n)
	      write(9,"(F10.4,2(3x,F12.8))") (n*binwidth)-binwidth/2.0, nrdf(sp1,sp2,n), integral
	    end do
	    close(9)
	  end do
	end do
	write(0,*) "Finished."
999	close(10)
	close(13)
	end program rdf

