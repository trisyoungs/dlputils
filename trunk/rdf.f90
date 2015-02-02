!	** rdf **
!	RDF program to calculate radial distribution functions between the centres-of-mass of
!	molecules. Automatically performed for all species pair combinations.
!	Calculates the density at each frame of the simulation (for NPT simulations)

	module rdfssdat
	  integer :: nbins
	  real*8 :: binwidth,boxvolume
	  real*8, allocatable :: purerdf(:),nrdf(:,:,:),sumhist(:,:,:)
	  integer, allocatable :: irdf(:,:,:)
	end module rdfssdat

	program rdf
	use dlprw; use rdfssdat; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	logical :: altheader = .FALSE., nonorm = .FALSE., zplus = .FALSE., zminus = .FALSE.
	integer :: n,m,a1,s1,s2,m1,m2,baselen,bin,nframes,success,o,nargs,numadded,framestodo = -1,frameskip = 0,framesdone, compairs(10,2), otherpairs(10,2)
	integer :: iargc
	real*8 :: dist,c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z,tx,ty,tz,bdens, integral

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
              n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") binwidth
              write(0,"(A,f6.2)") "Binwidth set to ",binwidth
            case ("-compair")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") s1
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(s1,1)
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(s1,2)
              write(0,"(A,3I4)") "Using COMpair for species ",s1, compairs(s1,:)
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
            case ("-nonorm")
              write(0,"(A)") "RDFs will not be normalised by species number density."
	      nonorm = .TRUE.
            case ("-otherpair")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") s1
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") otherpairs(s1,1)
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") otherpairs(s1,2)
              write(0,"(A,3I4)") "Using other pair (for surrounding molecules) for species ",s1, otherpairs(s1,:)
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
	
	call alloc_data(nspecies,nbins)

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
	do s1=1,nspecies
	  ! If compairs were specified, use that instead of COM
	  if (compairs(s1,1).ne.0) then
	    do m1=1,s_nmols(s1)
	      c1x = xpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,1)-1)
	      c1y = ypos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,1)-1)
	      c1z = zpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,1)-1)
	      c2x = xpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,2)-1)
	      c2y = ypos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,2)-1)
	      c2z = zpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,2)-1)
	      call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	      comx(s1,m1) = (c1x+tx)*0.5
	      comy(s1,m1) = (c1y+ty)*0.5
	      comz(s1,m1) = (c1z+tz)*0.5
	    end do
	  end if
	end do
	
	irdf = 0
	do s1=1,nspecies

	  do m1=1,s_nmols(s1)     ! Loop over all molecules of species 1...

	    ! Grab centre-of-mass (or compair) coordinate for central species
	    c1x=comx(s1,m1)
	    c1y=comy(s1,m1)
	    c1z=comz(s1,m1)
	
	    do s2=1,nspecies
	      numadded=0
	      ! Now loop over all molecules of second species....
	      do m2=1,s_nmols(s2)
		! Don't add if same species *and* same molecule
		if ((s1.eq.s2).and.(m1.eq.m2)) cycle
		! Grab the other centre coordinates if necessary
		if (otherpairs(s2,1).ne.0) then
	          c3x = xpos(s_start(s2)+(m2-1)*s_natoms(s2)+otherpairs(s2,1)-1)
	          c3y = ypos(s_start(s2)+(m2-1)*s_natoms(s2)+otherpairs(s2,1)-1)
	          c3z = zpos(s_start(s2)+(m2-1)*s_natoms(s2)+otherpairs(s2,1)-1)
	          c2x = xpos(s_start(s2)+(m2-1)*s_natoms(s2)+otherpairs(s2,2)-1)
	          c2y = ypos(s_start(s2)+(m2-1)*s_natoms(s2)+otherpairs(s2,2)-1)
	          c2z = zpos(s_start(s2)+(m2-1)*s_natoms(s2)+otherpairs(s2,2)-1)
	          call pbc(c2x,c2y,c2z,c3x,c3y,c3z,tx,ty,tz)
	          c2x = (c3x+tx)*0.5
	          c2y = (c3y+ty)*0.5
	          c2z = (c3z+tz)*0.5
		else
		  c2x=comx(s2,m2)
		  c2y=comy(s2,m2)
		  c2z=comz(s2,m2)
		end if
		! Get the shortest (MIM) distance between the two points
		call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
		! Check zplus or zminus
		if (zplus.and.(tz.gt.c1z)) cycle
		if (zminus.and.(tz.lt.c1z)) cycle
		dist=sqrt( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )
		! 'Add' this distance...
		bin=INT(dist/binwidth)+1
		if (bin.lt.nbins) irdf(s1,s2,bin) = irdf(s1,s2,bin) + 1
		numadded=numadded+1
	      end do
	      ! Consistency check - we should *always* add s_nmols() unless s1 = s2.
	      ! Don't do this check if zplus or zminus is in effect
	      if ((.not.zplus).and.(.not.zminus)) then
	        if (s1.eq.s2) then
	          if (numadded.NE.s_nmols(s1)-1) write(0,"(A,i6,a,i2,a,i2,' (',i7,')')") &
		    & "WARNING: Frame ",nframes," - Did not bin all molecules of ",s2," around ",s1,numadded
	        else
		  if (numadded.NE.s_nmols(s2)) write(0,*) &
		    & "WARNING: Did not bin all molecules...",nframes,s2,numadded
		end if
	      end if
	    end do
	  end do
	end do

	! Calc and store the normalised rdf....
	do s2=1,nspecies
	  if (.not.nonorm) call calcpurerdf(s2)
	  do s1=1,nspecies
	    if (nframes.EQ.1) then
	      if (s1.EQ.s2) write(0,"(A,I2,A,I2,A,I4,A)") "RDF of ",s2," about ",s1," : averaged over ",s_nmols(s1)-1," molecules."
	      if (s1.NE.s2) write(0,"(A,I2,A,I2,A,I4,A)") "RDF of ",s2," about ",s1," : averaged over ",s_nmols(s1)," molecules."
	    end if
	    do n=1,nbins
	      sumhist(s1,s2,n)=sumhist(s1,s2,n)+real(irdf(s1,s2,n))
	      if (s1.eq.s2) then
		if (nonorm) then
	          nrdf(s1,s2,n)=nrdf(s1,s2,n)+real(irdf(s1,s2,n)) / (s_nmols(s1)-1)
		else
	          nrdf(s1,s2,n)=nrdf(s1,s2,n)+real(irdf(s1,s2,n)) / (s_nmols(s1)-1) / purerdf(n)
		end if
	      else
		if (nonorm) then
		  nrdf(s1,s2,n)=nrdf(s1,s2,n)+real(irdf(s1,s2,n)) / s_nmols(s1)
		else
		  nrdf(s1,s2,n)=nrdf(s1,s2,n)+real(irdf(s1,s2,n)) / s_nmols(s1) / purerdf(n)
		end if
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

	do s1=1,nspecies
	  do s2=1,nspecies
	    call calcpurerdf(s2)
	    if (nonorm) then
	      if (zplus) then
		resfile=basename(1:baselen)//"u_zp_rdf"//CHAR(48+s1)//CHAR(48+s2)
	      else if (zminus) then
		resfile=basename(1:baselen)//"u_zm_rdf"//CHAR(48+s1)//CHAR(48+s2)
	      else
		resfile=basename(1:baselen)//"urdf"//CHAR(48+s1)//CHAR(48+s2)
	      end if
	    else
	      if (zplus) then
		resfile=basename(1:baselen)//"zp_rdf"//CHAR(48+s1)//CHAR(48+s2)
	      else if (zminus) then
		resfile=basename(1:baselen)//"zm_rdf"//CHAR(48+s1)//CHAR(48+s2)
	      else
		resfile=basename(1:baselen)//"rdf"//CHAR(48+s1)//CHAR(48+s2)
	      end if
	    end if
	    open(unit=9,file=resfile,form="formatted")
	    if (compairs(s1,1).eq.0) then
	      write(9,"(a,i2,a)") "# Species ",s1," calculated using all atoms for COM."
	    else
	      write(9,"(a,i2,a,i3,i3)") "# Species ",s1," calculated using two atoms for COM: ", compairs(s1,1), compairs(s1,2)
	    end if
	    if (compairs(s2,1).eq.0) then
	      write(9,"(a,i2,a)") "# Species ",s2," calculated using all atoms for COM."
	    else
	      write(9,"(a,i2,a,i3,i3)") "# Species ",s2," calculated using two atoms for COM: ", compairs(s2,1), compairs(s2,2)
	    end if
	    integral = 0.0
	    ! Normalise the RDFs with respect to the number of frames.
	    do n=1,nbins
	      integral = integral + sumhist(s1,s2,n) / framesdone / s_nmols(s1)
	      nrdf(s1,s2,n)=nrdf(s1,s2,n)/framesdone
	      ! write(9,"(F6.3,3x,F12.8)") (n*binwidth)-binwidth/2.0,nrdf(a1,n)
	      write(9,"(F10.4,2(3x,F12.8))") (n*binwidth)-binwidth/2.0,nrdf(s1,s2,n),integral
	    end do
	    close(9)
	  end do
	end do
	write(0,*) "Finished."
999	close(10)
	close(13)
	end program rdf

	subroutine calcpurerdf(sp)
	use dlprw; use rdfssdat
	implicit none
	! Calculate the 'pure' RDFs for supplied species
	integer :: sp,n,INFO
	real*8 :: const,numdens
	numdens = s_nmols(sp) / (cell(1)*cell(5)*cell(9))
	const = (4.0*3.141592654) / 3.0
	do n=1,nbins
	  purerdf(n) = (const*((n*binwidth)**3-((n-1)*binwidth)**3)) * numdens
	end do
	end subroutine calcpurerdf

	subroutine alloc_data(i,j)
	use rdfssdat; implicit none; integer :: n,i,j,status
	! i = nspecies, j = nbins
        allocate(irdf(i,i,j),stat=status); if (status.GT.0) stop "Allocation error for rdf()"
        allocate(purerdf(j),stat=status); if (status.GT.0) stop "Allocation error for purerdf()"
	allocate(nrdf(i,i,j),stat=status); if (status.GT.0) stop "Allocation error for nrdf()"
	allocate(sumhist(i,i,j),stat=status); if (status.GT.0) stop "Allocation error for sumhist()"
	! Initialise the arrays...
	nrdf = 0.0
	sumhist = 0.0
	irdf = 0
	purerdf = 0.0
	end subroutine alloc_data
