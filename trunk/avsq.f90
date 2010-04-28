!	** avsq **
!	Average the supplied files into a single dataset

	program avsq
	use dlprw; use utility; use parse
	implicit none

	real*8, parameter :: pi = 3.14159265358979
	character*80, allocatable :: inputfiles(:)
	character*80 :: outfile
	real*8, allocatable :: sq(:,:), newsq(:), error(:)
	real*8 :: binwidth
	integer :: n, b, nfiles, nbins, nargs, iargc
	logical :: success

	nargs = iargc()
	if (nargs.lt.2) then
	  write(0,*) "Usage : avsq <files.....> <outputfile>"
	  stop
	end if
	nfiles = nargs-1

	! Allocate filename array
	allocate(inputfiles(nfiles))
	write(0,*) nfiles, "files supplied..."

	do n=1,nfiles
	  call getarg(n,inputfiles(n))
	end do
	call getarg(nargs,outfile)

	! Open output file
	open(unit=15,file=outfile, form="formatted", status='new')
 
	! Open first file and determine number of bins
	open(unit=11,file=inputfiles(1), form="formatted", status="old")
	! First few lines are comments (beginning with '#', and giving 0 arguments)
50	success = readline(11)
	if (nargsparsed.eq.0) goto 50
	binwidth = argr(1) * 2.0
	write(0,*) "Bin width in file is ", binwidth
	nbins=1
101	success = readline(11)
	if (nargsparsed.lt.1) goto 110
	nbins = nbins + 1
	goto 101

110	close(11)
	write(0,*) "Number of bins in file is ", nbins

	! Allocate data arrays
	allocate(sq(nfiles,nbins))
	allocate(newsq(nbins))
	allocate(error(nbins))

	! Read in all data
	do n=1,nfiles
	  open(unit=11,file=inputfiles(n), form="formatted", status="old")
	  write(0,*) "Reading file :", inputfiles(n)
	  success = readline(11)
	  success = readline(11)
	  do b=1,nbins
	    success = readline(11)
	    ! Check for expected binwidth for first point
	    if (b.eq.1) then
	      if (abs(binwidth-argr(1).gt.0.0001)) stop "Incompatible binwidths."
	    end if
	    sq(n,b) = argr(2)
	  end do
	  close(11)
	end do

	! Calculate average and S.D. error
	newsq = 0.0
	do n=1,nfiles
	  newsq = newsq + sq(n,:)
	end do
	newsq = newsq / nfiles
	
	error = 0.0
	do n=1,nfiles
	  do b=1,nbins
	    error(b) = error(b) + (sq(n,b)-newsq(b))**2
	  end do
	end do
	error = error / nfiles
	error = sqrt(error)

	! Write new data
	do b=1,nbins
	  write(15,"(3f12.6)") (b-0.5)*binwidth, newsq(b), error(b)
	end do
	close(15)

	end program avsq
