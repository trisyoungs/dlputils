	! Program to sum the results in 2 or more rdf output files, and average them to create
	! a new file.
	
	program rdfsum
	use parse
	implicit none
	character*80, allocatable :: files(:)
	character*80 :: discard
	integer :: nargs,sndata,n,m,count
	integer :: iargc
	real*8, allocatable :: rdf(:,:)  ! Will be (N,3) : (,1) = bin centre, (,2) = rdf value, (,3) = integral
	real*8 :: a,c,i,x

        nargs = iargc()
        if (nargs.EQ.0) stop " Usage : rdfsum <file1> <file2> ... <fileN>  (results written to stdout)"
	if (nargs.EQ.1) write(0,*) "Warning! Only one file supplied - results will be a copy."
	allocate(files(nargs))
        do n=1,nargs
          call getarg(n,files(n))
        end do

!	From the first of these files, determine the number of data points we need to read in...
	open(unit=15,file=files(1),form='formatted',status='old')
	count = 0
	if (discard(1:4).NE." Bin") rewind(15)
	! Read through the data until we reach the end of the file...
10	if (.not.readline(15)) goto 20
	if (nargsparsed.ne.0) count = count + 1
	goto 10
20	close(15)

	! Can now assign the array
	allocate(rdf(count,3)) 
	rdf = 0.0

	! Process the files given
	do n=1,nargs
	  open(unit=15,file=files(n),form='formatted',status='old')
	  ! Header check...
	  read(15,"(A80)") discard
	  if (discard(1:4).NE." Bin") rewind(15)
	  do m=1,count
	    if (.not.readline(15)) exit
! 	    read(15,"(F6.3,3x,F12.8)",end=20,err=20) a,c,i
	    rdf(m,1) = argr(1)
	    rdf(m,2) = rdf(m,2) + argr(2)
	    rdf(m,3) = rdf(m,3) + argr(3)
	  end do
	  close(15)
	end do

	! Take the average of rdf(:,2)
	rdf(:,2) = rdf(:,2) / real(nargs)
	write(0,"(A,I2,A)") "Averaged RDF of ",nargs," files"

	! Print out the data
	do n=1,count
	  write(6,"(F6.3,3x,F12.8,3x,e12.6)") rdf(n,1),rdf(n,2),rdf(n,3)
	end do

	end program rdfsum
