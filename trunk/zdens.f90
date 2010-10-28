!	** zdens **
!	Calculate density profile along z-direction of unit cell

	program zdens
	use dlprw; use utility
	implicit none
	integer :: nbins
	real*8 :: binwidth, boxvolume, z, zbase
	real*8, allocatable :: dens(:,:)
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	logical :: altheader = .FALSE.
	integer :: n,m,a1,s1,m1,baselen,bin,nframes,success,nargs,framestodo = -1, framestoskip=-1
	integer :: iargc, nskipped
	integer, allocatable :: comatom(:)

	binwidth=0.1   ! In Angstroms
	nargs = iargc()
	if (nargs.lt.3) stop "Usage : zdens <DLP HISTORYfile> <DLP OUTPUTfile> <base z> [-bin width] [-header hisfile] [-comatom sp atom] [-skip n] [-frames n]"
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
            case ("-skip")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestoskip
              write(0,"(A,I6)") "Frames to skip: ",framestoskip
            case ("-comatom")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") s1
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") comatom(s1)
              write(0,"(A,I2,A,I2,A)") "Will use atom ",comatom(s1)," for species ",s1," COM"
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
	
	allocate(dens(nspecies,-nbins-1:nbins+1))
	dens = 0.0

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
	nskipped=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.EQ.-1) goto 799  ! File error....
	nskipped = nskipped + 1
	if (nskipped.le.framestoskip) goto 101
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

	do s1=1,nspecies
	  resfile=basename(1:baselen)//"zdn"//CHAR(48+s1)
	  open(unit=9,file=resfile,form="formatted",status="replace")
	  write(9,"(a,f15.4)") "# Z relative to ",zbase
	  ! Normalise the RDFs with respect to the number of frames
	  do n=-nbins-1,nbins+1
	    dens(s1,n) = dens(s1,n) / nframes
	    write(9,*) n * binwidth, dens(s1,n)
	  end do
	  close(9)
	end do
999	write(0,*) "Finished."
	end program zdens