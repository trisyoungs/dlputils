!	** zdens **
!	Calculate density profile along z-direction of unit cell

	program zdens
	use dlprw; use utility
	implicit none
	integer :: nbins
	real*8 :: binwidth, boxvolume, z, zbase, centrez
	real*8, allocatable :: dens(:,:), avgdens(:,:)
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	logical :: altheader = .FALSE., average = .FALSE.
	integer :: n,m,a1,s1,m1,baselen,bin,nframes,success,nargs,framestodo = -1, framestodiscard=-1
	integer :: iargc, nskipped, centrebin, lowbin, hibin
	integer, allocatable :: comatom(:)

	binwidth=0.1   ! In Angstroms
	nargs = iargc()
	if (nargs.lt.3) stop "Usage : zdens <HISTORYfile> <OUTPUTfile> <base z> [-bin width] [-header hisfile] [-comatom sp atom] [-discard n] [-frames n] [-symm cz]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(f15.4)") zbase
	
	! Open and check the files...
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
	allocate(comatom(nspecies))
	comatom = 0

	n = 2
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-bin")
              n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") binwidth
              write(0,"(A,f10.4)") "Binwidth set to ",binwidth
            case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
	      altheader = .TRUE.
            case ("-frames")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
              write(0,"(A,I6)") "Frames to do: ",framestodo
            case ("-discard")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodiscard
              write(0,"(A,I6)") "Frames to discard: ",framestodiscard
            case ("-comatom")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") s1
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") comatom(s1)
              write(0,"(A,I2,A,I2,A)") "Will use atom ",comatom(s1)," for species ",s1," COM"
            case ("-symm")
              n = n + 1; call getarg(n,temp); read(temp,"(f20.14)") centrez
	      average = .TRUE.
              write(0,"(A,f10.5)") "Distributions will be symmetrized about z = ", centrez
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

	nbins = maxval(cell) / binwidth
	if (nargs.EQ.3) then
	  write(0,*) "Using default binwidth of 0.1 Angstroms"
	else
	  write(0,"(A,F6.3,A)") "Using specified binwidth of ",binwidth," Angstroms"
	end if
	write(0,"(A,I5,A)") "There will be ",nbins," histogram bins."

	if (average) then
	  avgdens = 0.0
	  ! Determine central bin and check its in range
	  centrebin = centrez/binwidth
	  if ((centrebin.lt.(-nbins-1)).or.(centrebin.gt.(nbins+1))) stop "Supplied centreZ is out of bin range."
	end if
	
	allocate(dens(nspecies,-nbins-1:nbins+1))
	allocate(avgdens(nspecies,-nbins-1:nbins+1))
	dens = 0.0

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
	nskipped=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nskipped = nskipped + 1
	if (nskipped.le.framestodiscard) goto 101
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	call calc_com

	m = 0
	do s1=1,nspecies
	  do m1=1,s_nmols(s1)     ! Loop over all molecules of species
	    ! Select COM Z coordinate
	    if (comatom(s1).eq.0) then
	      z = comz(s1,m1)
	    else
	      z = zpos(m+comatom(s1))
	    endif
	    ! Normalise z coordinate relative to 'zbase'
	    ! !!!Fold in z direction if greater than cell(9)*0.5
	    z = z - zbase
	    if (z.lt.0.0) z = z + cell(9)
	    !if (z.lt.(-cell(9)*0.5)) z = z + cell(9)
	    !if (z.gt.(cell(9)*0.5)) z = z - cell(9)
	    bin = (z / binwidth) + 1
	    dens(s1,bin) = dens(s1,bin) + 1
	    ! Increase atom counter
	    m = m + s_natoms(s1)
	  end do
	end do

	if (nframes.eq.framestodo) goto 801
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

	write(0,*) "Read ",nframes,"frames."
	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="density"
	  baselen=7
	else
	  basename=hisfile(1:baselen)
	endif

	! Average data about centrepoint of distribution
	if (average) then
	  avgdens = 0.0
	  do s1=1,nspecies
	    avgdens(s1,centrebin) = dens(s1,centrebin)
	    do n=1,nbins
	      lowbin = centrebin-n
	      if (lowbin.lt.(-nbins-1)) lowbin = lowbin + (2*nbins+3)
	      hibin = centrebin+n
	      if (hibin.gt.(nbins+1)) hibin = hibin - (2*nbins+3)
	      avgdens(s1,lowbin) = 0.5*(dens(s1,lowbin)+dens(s1,hibin))
	      avgdens(s1,hibin) = avgdens(s1,lowbin)
	    end do
	  end do
	end if

	! Normalise against number of frames
	dens = dens / nframes
	avgdens = avgdens / nframes

	! Write data
	do s1=1,nspecies
	  resfile=basename(1:baselen)//"zdn"//CHAR(48+s1)
	  open(unit=9,file=resfile,form="formatted",status="replace")
	  write(9,"(a,f15.4)") "# Z relative to ",zbase
	  do n=-nbins-1,nbins+1
	    write(9,"(3e15.5)") (n+0.5) * binwidth, dens(s1,n), avgdens(s1,n)
	  end do
	  close(9)
	end do
999	write(0,*) "Finished."
	stop
	end program zdens
