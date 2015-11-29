	! Apply a gaussian filter to a pdens file

	program pdensgauss
	use parse; use pdensrw
	implicit none
	character*80 :: infile, outfile
	character*20 :: looporder, temp
	character :: c
	real*8, parameter :: pi = 3.14159265358979
	integer :: nargs, n, x, y, z, point(3), extents(3), i, j, k
	logical :: success, minLimit = .false., maxLimit = .false.
	real*8 :: v(3), sigma, mag(3), xyz, twosigma2, minValue, maxValue
	real*8, allocatable :: gaussvalues(:,:,:)
	type(PDens) :: original, gauss
	integer :: iargc

        nargs = iargc()
        if (nargs.lt.3) stop "Usage: pdensgauss <pdensfile> <outputfile> <sigma> [-min value] [-max value]"
        call getarg(1,infile)
        call getarg(2,outfile)
	sigma = getargr(3)

	! Set defaults
	n = 3
        do
          n = n + 1; if (n.gt.nargs) exit
          call getarg(n,temp)
          select case (temp)
	    case ("-min")
	      n = n + 1; minValue = getargr(n)
	      minLimit = .true.
	      write(0,*) "Input grid will be pruned, removing values below ", minValue
	    case ("-max")
	      n = n + 1; maxValue = getargr(n)
	      maxLimit = .true.
	      write(0,*) "Input grid will be pruned, removing values above ", maxValue
	    case default
	      write(0,*) "Unrecognised argument :",temp
	      stop
	  end select
	end do

	! Open original pdens file
	if (.not.loadPDens(infile,original)) stop "Couldn't load original pdens file."

	! Copy necessary quantities over to new pdens, and allocate its array
	gauss%loop = original%loop
	gauss%origin = original%origin
	gauss%axes = original%axes
	if (.not.allocPDens(gauss, original%gridMin(1),original%gridMin(2),original%gridMin(3),  &
 		& original%gridMax(1),original%gridMax(2),original%gridMax(3))) stop "Failed to allocate new PDens."

	! Prune input grid
	if (minLimit.or.maxLimit) then
	  write(0,*) "Pruning input grid..."
	  ! Loop over original data
	  do x=original%gridMin(1),original%gridMax(1)
	    write(0,"(a,i4,a,i4)") "At x = ", x, ", max = ", original%gridMax(1)
	    do y=original%gridMin(2),original%gridMax(2)
	      do z=original%gridMin(3),original%gridMax(3)
		if (minLimit.and.(original%grid(x,y,z).lt.minValue)) then
		  original%grid(x,y,z) = 0.0
		  cycle
		end if
		if (maxLimit.and.(original%grid(x,y,z).gt.maxValue)) then
		  original%grid(x,y,z) = 0.0
		  cycle
		end if
	      end do
	    end do
	  end do
	end if

	! Determine number of gridpoints to step in each direction, based on gaussian width
	extents = 0
	do n=1,3
	  mag(n) = sqrt(original%axes((n-1)*3+1)*original%axes((n-1)*3+1) + original%axes((n-1)*3+2)*original%axes((n-1)*3+2) + original%axes((n-1)*3+3)*original%axes((n-1)*3+3))
	  do while ( exp(  (-((mag(n)*extents(n))*(mag(n)*extents(n))))  /  (2*sigma*sigma)) .gt.1.0e-2 )
	    write(0,"(a,i1,a,i2,a,f10.4)") "Axis ", n, " extent of ", extents(n), " gives Gaussian value of ", exp(  (-((mag(n)*extents(n))*(mag(n)*extents(n))))  /  (2*sigma*sigma))
	    extents(n) = extents(n) + 1
	  end do
	end do
	extents = extents - 1
	write(0,*) "Extents for gaussian filtering (in gridpoints): ", extents
	if (sum(extents).eq.0) stop "Extents are zero (nothing to do). Try increasing sigma."

	! Work out list of surrounding points to consider at each gridpoint
	twosigma2 = 2.0 * sigma * sigma
	allocate(gaussvalues(-extents(1):extents(1),-extents(2):extents(2),-extents(3):extents(3)))
	do i=-extents(1),extents(1)
	  do j=-extents(2),extents(2)
	    do k=-extents(3),extents(3)
	      xyz = (mag(1)*i)*(mag(1)*i) + (mag(2)*j)*(mag(2)*j) + (mag(3)*k)*(mag(3)*k)
	      gaussvalues(i,j,k) = (1.0 / (4.0*pi*pi*sigma*sigma*sigma*sigma)) * exp(-xyz/twosigma2)
	    end do
	  end do
	end do
	write(0,*) gaussvalues
	
	! Loop over original data
	do x=original%gridMin(1),original%gridMax(1)
	  write(0,"(a,i4,a,i4)") "At x = ", x, ", max = ", original%gridMax(1)
	  do y=original%gridMin(2),original%gridMax(2)
	    do z=original%gridMin(3),original%gridMax(3)

	      ! Loop over gaussian extent
	      do i=-extents(1),extents(1)
	        do j=-extents(2),extents(2)
	          do k=-extents(3),extents(3)
		    point(1) = x+i
		    if ((point(1).lt.original%gridMin(1)).or.(point(1).gt.original%gridMax(1))) cycle
		    point(2) = y+j
		    if ((point(2).lt.original%gridMin(2)).or.(point(2).gt.original%gridMax(2))) cycle
		    point(3) = z+k
		    if ((point(3).lt.original%gridMin(3)).or.(point(3).gt.original%gridMax(3))) cycle

		    !xyz = (mag(1)*i)*(mag(1)*i) + (mag(2)*j)*(mag(2)*j) * (mag(3)*k)*(mag(3)*k)
		    !g = grid(x,y,z) * exp(-xyz/twosigma2)
		    gauss%grid(point(1),point(2),point(3)) = gauss%grid(point(1),point(2),point(3)) + original%grid(x,y,z)*gaussvalues(i,j,k)
		  end do
		end do
	      end do

	    end do
	  end do
	end do

	! Normalise
	!gaussgrid = gauss%grid / (4.0*pi*pi*sigma*sigma*sigma*sigma)
	gauss%grid = gauss%grid * (sum(original%grid) / sum(gauss%grid))
	!gauss%grid = gauss%grid * (maxval(original%grid) / maxval(gauss%grid))
	write(0,"(a,2e12.5)") "Sanity sums are ", sum(original%grid) , sum(gauss%grid)
	
	! Save data
	if (.not.savePDens(outfile, gauss)) stop "Failed to save smoothed data."

	end program pdensgauss
