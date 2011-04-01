	! Program to sum the results in 2 or more acf output files, and average them to create the total
	
	program acfsum
	use parse
	implicit none
	character*80, allocatable :: files(:)
	character*80 :: discard
	integer :: nargs,n,m,count,i
	logical :: success
	integer :: iargc
	real*8, allocatable :: acf(:,:), t(:), accum(:)
	real*8 :: acc

        nargs = iargc()
        if (nargs.EQ.0) stop " Usage : acfsum <file1> <file2> ... <fileN>  (results written to stdout)"
	if (nargs.EQ.1) write(0,*) "Warning! Only one file supplied - results will be a copy."
	allocate(files(nargs))
        do n=1,nargs
          call getarg(n,files(n))
        end do

!	From the first of these files, determine the number of data points we need to read in...
	open(unit=15,file=files(1),form='formatted',status='old')
	! Skip first line since it is the commented column header
	success = readline(15)
	! Read through the data until we reach the end of the file...
	count = 0
10	if (.not.readline(15)) goto 20
	if (nargsparsed.ne.0) count = count + 1
	goto 10
20	close(15)

	! Can now assign the arrays
	allocate(acf(count,7)) 
	allocate(accum(count))
	allocate(t(count))

	acf = 0.0
	accum = 0.0
	t = 0.0

	! Process the files given
	do n=1,nargs
	  open(unit=15,file=files(n),form='formatted',status='old')
	  ! Discard header
	  success = readline(15)
	  do m=1,count
	    if (.not.readline(15)) exit
! 	    read(15,"(F6.3,3x,F12.8)",end=20,err=20) a,c,i
	    t(m) = argr(1)
	    acc = argr(9)
	    accum(m) = accum(m) + acc
	    do i=1,7
	      acf(m,i) = acf(m,i) + argr(i+1) * acc
	    end do
	  end do
	  close(15)
	end do

	! Form averages
	do m=1,count
	  acf(m,1:7) = acf(m,1:7) / accum(m)
	end do
	write(0,"(A,I2,A)") "Averaged ACF over ",nargs," files"

	! Print out the data
	write(6,"(10a14)") "#t","total","xx","yy","zz","xy","xz","yz","acc"
	do m=1,count
	  write(6,"(9f14.6,e14.6)") t(m), (acf(m,n),n=1,7), accum(m)
	end do

	end program acfsum
