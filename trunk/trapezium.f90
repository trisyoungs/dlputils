	program trapezium
	use parse
	implicit none
	
	real*8, allocatable :: a(:), b(:), integral(:)
	integer :: ndata = 0, ncols = 0, n
 
	! Read line from standard input
100	if (.not.readline(5)) goto 200
	if (nargsparsed.eq.0) goto 100
	if (nargsparsed.eq.-1) goto 200
	ndata = ndata + 1
	if (ndata.eq.1) then
	  allocate(a(nargsparsed))
	  allocate(b(nargsparsed))
	  allocate(integral(nargsparsed))
	  ncols = nargsparsed
	end if
	! Get data from args()
	do n=1,ncols
	  a(n) = argr(n)
	end do
	! Assume a(1) is the abcissa
	if (ndata.gt.1) then
	  do n=2,ncols
	    integral(n) = integral(n) + ( a(1)-b(1) )*( a(n)+b(n) )*0.5
	  end do
	end if
	! Store this data as the old data (b)
	b = a
	write(6,"('x=',11e15.5)") a(1), integral(2:ncols)
	goto 100


200	write(6,*) "Read ",ndata,"points per column from file"
	write(6,*) "Final integrals:"
	write(6,"('    ',10e15.5)")  integral(2:ncols)

	end program trapezium
