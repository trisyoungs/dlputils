!	** acfcat **
!	Concatenate quantity files

	program acfcat
	use dlprw; use utility
	implicit none

	character*80 :: files(100),outfile,arg
	character*20 :: temp
	integer :: n,nfiles, nframes
	integer :: iargc, qmax, nargs
	real*8, allocatable :: qx(:), qy(:), qz(:)

	nargs = iargc()
	if (nargs.lt.2) then
	  write(0,"(a)") "Usage : acfcat <nmols> <file1> <file2> ... <fileN> <target>"
	  stop
	end if
	call getarg(1,temp); read(temp,"(i6)") qmax
	write(0,*) "Number of molecules in file (qmax) = ", qmax
	nfiles = 0
	do n=2,nargs
	  nfiles = nfiles + 1
	  call getarg(n, files(nfiles))
	  write(6,"(a,i3,a,a80)") "File ", nfiles, " is : ", files(nfiles)
	end do

	! Last specified filename is output filename
	outfile = files(nfiles)
	nfiles = nfiles - 1
	if (nfiles.le.0) stop "Error: No input files were given." 
	write(6,*) "Output file will be : ", outfile
	! Test open the output file....
	open(unit=20, file=outfile, status='old', err=5)
	stop "Will not overwrite existing output file."
	
5	close(20)
	write(6,*) "Concatenating..."

	allocate (qx(qmax))
	allocate (qy(qmax))
	allocate (qz(qmax))

	! Open output file (binary)
	open(unit=20, file=outfile, status='new', form='unformatted')

	nframes=0

	do n=1,nfiles
	  ! Open nth file
	  write(0,*) files(n)
	  open(unit=11, file=files(n), status='old', form='unformatted', err=799)
	  ! Read data for all molecules
100	  read(11,err=798,end=200) qx
	  read(11,err=798,end=797) qy
	  read(11,err=798,end=797) qz
	  ! Write data for all molecules
	  write(20) qx
	  write(20) qy
	  write(20) qz
	  nframes = nframes + 1
	  if (mod(nframes,100).eq.0) write(0,*) nframes
	  goto 100
200	  close(11)
	end do

	goto 801
797	write(0,*) "End of file while reading qy or qz from:", files(nfiles)
	goto 801
798	write(0,*) "Error reading quantity from:", files(nfiles)
	goto 801
799	write(0,*) "Error opening file:", files(nfiles)
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",nframes,")"
	goto 801
801	write(0,*) ""

	close(15)

999	write(0,*) "Finished."
	end program acfcat
