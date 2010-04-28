	program args
	implicit none
	character*80 :: arg, infile1, infile2, outfile
	integer :: iargc, nargs, n

	nargs = iargc()

	write(6,*) nargs

	do n=1,nargs

	  call getarg(n,arg)
	  write(6,*) "Argument", n, " is ",arg

	end do

	call getarg(1,infile1)
	call getarg(2,infile2)
	call getarg(3,outfile)


	open(unit=10, file=infile1, status="old")

	end program args

