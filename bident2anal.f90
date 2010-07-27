!	** bident2anal **
!	Analyse contacts placed in bident2 dump file

	module dat
	integer :: ncentre, nouter, ncentresites
	integer, allocatable :: outerdata(:,:,:,:), centredata(:,:,:,:,:), ncentredata(:,:,:), nouterdata(:,:)
	integer, allocatable :: maskarraysame(:,:), maskarrayother(:,:)
	end module dat

	program bident2anal
	use parse; use dat
	implicit none
	character*80 :: dumpfile, basename
	character*20 :: temparg
	integer, parameter :: MAXSITES = 20, MAXCONTACTS = 20
	integer :: nframes,nargs,n,m1,m2,site,baselen,framestodo,f,m,c
	integer :: iargc, bitmask
	logical :: success

	nargs = iargc()
	if (nargs.lt.2) then
	  write(0,"(a120)") "Usage : bident2anal <dumpfile> <framestodo> [selectors...]"
	  stop
	end if
	call getarg(1,dumpfile)
	call getarg(2,temparg); read(temparg,"(i10)") framestodo
	
	! Open and check the files...
	open(unit=11,file=dumpfile,form='formatted',status='old')

	n = 2
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temparg)
          select case (temparg)
            case ("-frames")
              !n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") framestodo
	    case default
	      write(0,*) "Unrecognised argument:", temparg
	      stop
	  end select
	end do
	
	! Read in preliminary data to get array sizes
	success = readline(11); ncentre = argi(2)
	success = readline(11); nouter = argi(2)
	success = readline(11); ncentresites = argi(2)
	
	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (dumpfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="rdfresults."
	  baselen=11
	else
	  basename=dumpfile(1:baselen)
	endif

	write(0,*) "ncentre = ", ncentre
	write(0,*) "nouter = ", nouter
	write(0,*) "ncentresites = ", ncentresites
	
	allocate(centredata(framestodo, ncentre, ncentresites, MAXCONTACTS, 0:2))
	allocate(ncentredata(framestodo, ncentre, ncentresites))
	allocate(outerdata(framestodo, nouter, MAXCONTACTS, 0:3))
	allocate(nouterdata(framestodo, nouter))
	allocate(maskarraysame(ncentresites,2**(ncentresites*2)))
	allocate(maskarrayother(ncentresites,2**(ncentresites*2)))

	centredata = 0
	ncentredata = 0
	outerdata = 0
	nouterdata = 0

	! Read in data from file
100	nframes=0
101	success = readline(11)
	if (.not.success) goto 801  ! End of file encountered....
	select case (arg(1))
	  case ("frame")
	    nframes = nframes + 1
	    if (nframes.gt.framestodo) goto 800
	    if (mod(nframes,100).EQ.0) write(0,*) nframes
	  case ("centre")
	    m1 = argi(2)
	    site = argi(4)
	    m2 = argi(6)
	    bitmask = argi(8)
	    n = ncentredata(nframes, m1, site) + 1
	    ncentredata(nframes, m1, site) = n
	    centredata(nframes, m1, site, n, 1) = m2
	    centredata(nframes, m1, site, n, 2) = bitmask
	    n = nouterdata(nframes, m2) + 1
	    nouterdata(nframes, m2)  = n
	    outerdata(nframes, m2, n, 1) = m1
	    outerdata(nframes, m2, n, 2) = site
	    outerdata(nframes, m2, n, 3) = bitmask
	  case ("outer")
	    !m1 = argi(2)
	    !n = nouterdata(nframes, m1) + 1
	    !nouterdata(nframes, m1) = n
	    !outerdata(nframes, m1, n, 1) = argi(4)
	    !outerdata(nframes, m1, n, 2) = argi(6)
	end select

	! Next data
	goto 101

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999
798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended before framestodo was fulfilled..."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",framestodo,")"
	goto 801
800	write(0,*) "Framestodo was fulfilled."
801	write(0,*) ""

	close(11)

	! 1) For bidentate outer species, what is the contact pattern to other molecules
	call reset()
	do f=1,nframes
	  do m2 = 1,5
	    do n=1,nouterdata(f,m2)
	      m1 = outerdata(f,m2,n,1)
	      site = outerdata(f,m2,n,2)
	      ! Is this interaction bidentate?
	      if (outerdata(f,m2,n,3).eq.3) exit
	    end do
	    ! If we found one, construct a bitmask
	    if (n.le.nouterdata(f,m2)) then
	      ! Construct bitmask of other interactions
	      bitmask = 0
	      c = 0
	      do m=1,nouterdata(f,m2)
		if (n.eq.m) cycle
		bitmask = bitmask + 2**(outerdata(f,m2,m,3)-1+c)
	write(0,*) c,outerdata(f,m2,m,3), 2**(outerdata(f,m2,m,3)-1+c)
		c = c + 3
	      end do
	write(0,*) "Final bitmask", bitmask
	    end if
	  end do
	end do

	write(0,*) "Finished."
999	close(14)
	end program bident2anal

	subroutine reset()
	use dat
	implicit none
	maskarraysame = 0
	maskarrayother = 0
	centredata(:,:,:,:,0) = 0
	outerdata(:,:,:,0) = 0
	end subroutine reset
