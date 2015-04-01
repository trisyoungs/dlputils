	! Apply a gaussian filter to a pdens file

	program pdensgauss
	use parse; use pdensrw
	implicit none
	character*80 :: infile, outfile
	character*20 :: looporder, temp
	character :: c
	real*8, parameter :: pi = 3.14159265358979
	integer :: nargs, n, x, y, z, point(3), extents(3), i, j, k
	logical :: gnufile = .false., success
	real*8 :: v(3), sigma, mag(3), xyz, twosigma2, g
	real*8, allocatable :: gaussvalues(:,:,:)
	type(PDens) :: original, gauss
	integer :: iargc

        nargs = iargc()
        if (nargs.lt.3) stop "Usage: pdensgauss <pdensfile> <outputfile> <sigma>"
        call getarg(1,infile)
        call getarg(2,outfile)
	call getarg(3,temp); read(temp,"(f20.14)") sigma

	! Set defaults
	n = 3
        do
          n = n + 1; if (n.gt.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-gnuplot")
	      gnufile = .true.
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
	if (.not.allocPDens(gauss, original%ngrid)) stop "Failed to allocate new PDens."

	! Determine number of gridpoints to step in each direction, based on gaussian width
	extents = 1
	do n=1,3
	  mag(n) = sqrt(original%axes((n-1)*3+1)*original%axes((n-1)*3+1) + original%axes((n-1)*3+2)*original%axes((n-1)*3+2) + original%axes((n-1)*3+3)*original%axes((n-1)*3+3))
	  do while ( exp(  (-((mag(n)*extents(n))*(mag(n)*extents(n))))  /  (2*sigma*sigma)) .gt.1.0e-2 )
	    extents(n) = extents(n) + 1
	  end do
	end do
	write(0,*) "Extents for gaussian filtering (in gridpoints): ", extents

	! Work out list of surrounding points to consider at each gridpoint
	twosigma2 = 2.0 * sigma * sigma
	allocate(gaussvalues(-extents(1):extents(1),-extents(2):extents(2),-extents(3):extents(3)))
	do i=-extents(1),extents(1)
	  do j=-extents(2),extents(2)
	    do k=-extents(3),extents(3)
	      xyz = (mag(1)*i)*(mag(1)*i) + (mag(2)*j)*(mag(2)*j) * (mag(3)*k)*(mag(3)*k)
	      gaussvalues(i,j,k) = exp(-xyz/twosigma2)
	    end do
	  end do
	end do
	
	! Loop over original data
	do x=-original%ngrid,original%ngrid
	  write(0,"(a,i4,a,i4)") "At x = ", x+original%ngrid+1, " of ", 2*original%ngrid+1
	  do y=-original%ngrid,original%ngrid
	    do z=-original%ngrid,original%ngrid

	      ! Loop over gaussian extent
	      do i=-extents(1),extents(1)
	        do j=-extents(2),extents(2)
	          do k=-extents(3),extents(3)
		    point(1) = x+i
		    if ((point(1).lt.-original%ngrid).or.(point(1).gt.original%ngrid)) cycle
		    point(2) = y+j
		    if ((point(2).lt.-original%ngrid).or.(point(2).gt.original%ngrid)) cycle
		    point(3) = z+k
		    if ((point(3).lt.-original%ngrid).or.(point(3).gt.original%ngrid)) cycle

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
	!gaussgrid = gaussgrid / (4.0*pi*pi*sigma*sigma*sigma*sigma)
	!gaussgrid = gaussgrid * (sum(grid) / sum(gaussgrid))
	gauss%grid = gauss%grid * (maxval(original%grid) / maxval(gauss%grid))
	
	! Save data
	if (.not.savePDens(outfile, gauss)) stop "Failed to save smoothed data."

	end program pdensgauss
