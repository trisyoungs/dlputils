	! Program to sum the results in 2 or more acf output files, and average them to create the total
	
	program acfsum
	use parse
	implicit none
	character*80, allocatable :: files(:)
	character*80 :: discard
	integer :: nargs,n,i, ndata, nexpected
	logical :: success
	integer :: iargc
	real*8, allocatable :: acf(:,:), t(:), accum(:)
	real*8 :: acc, val(7), time

        nargs = iargc()
        if (nargs.EQ.0) stop " Usage : acfsum <file1> <file2> ... <fileN>  (results written to stdout)"
	if (nargs.EQ.1) write(0,*) "Warning! Only one file supplied - results will be a copy."
	
	write(0,"(a,i5,a)") "Performing average of ", nargs, " files"
	allocate(files(nargs))
        do n=1,nargs
          call getarg(n,files(n))
        end do

!	From the first of these files, determine the number of data points we need to read in...
	open(unit=15,file=files(1),form='formatted',status='old')
	! Skip first line since it is the commented column header
	success = readline(15)
	! Read through the data until we reach the end of the file...
	nexpected = 0
10	if (.not.readline(15)) goto 20
	if (nargsparsed.ne.0) nexpected = nexpected + 1
	goto 10
20	close(15)
	write(0,*) "Number of points in the first file : ", nexpected

	! Can now assign the arrays
	allocate(acf(nexpected,7)) 
	allocate(accum(nexpected))
	allocate(t(nexpected))

	acf = 0.0
	accum = 0.0
	t = 0.0

	! Process the files given
	do n=1,nargs
	  open(unit=15,file=files(n),form='formatted',status='old')
	write(0,*) "Processing file ",n
	  ! Discard header
	  success = readline(15)
	  ndata = 0
	  !if (.not.readline(15)) exit
25	  read(15,*,end=50,err=50) time, val(1), val(2), val(3), val(4), val(5), val(6), val(7), acc
	  ndata = ndata + 1
	  if (ndata.gt.nexpected) then
	    write(0,*) "WARNING - This file contains more data than the first...\n"
	    goto 50
	  end if
	  t(ndata) = time
	  accum(ndata) = accum(ndata) + acc
	  do i=1,7
	    acf(ndata,i) = acf(ndata,i) + val(i) * acc
	  end do
	  goto 25
50	  close(15)
	end do

	! Form averages
	do n=1,nexpected
	  acf(n,1:7) = acf(n,1:7) / accum(n)
	end do
	write(0,"(A,i3,A)") "Averaged ACF over ",nargs," files"

	! Print out the data
	write(6,"(10a14)") "#t","total","xx","yy","zz","xy","xz","yz","acc"
	do n=1,nexpected
	  write(6,"(9f14.6,e14.6)") t(n), (acf(n,i),i=1,7), accum(n)
	end do

	end program acfsum
