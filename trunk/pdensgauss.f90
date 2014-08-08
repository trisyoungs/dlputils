	! Apply a gaussian filter to a pdens file

	program pdensgauss
	use parse
	implicit none
	character*80 :: infile, outfile
	character*20 :: looporder, temp
	character :: c
	real*8, parameter :: pi = 3.14159265358979
	integer :: ngrid(3), nargs, n, x, y, z, loop(3), point(3), extents(3), i, j, k
	logical :: gnufile = .false., success
	real*8 :: axes(9), origin(3), v(3), sigma, mag(3), xyz, twosigma2, g
	real*8, allocatable :: grid(:,:,:), gaussgrid(:,:,:)
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

	! Open file
	open(unit=11, file=infile, form='formatted', status='old')

	! Format of data is :
	! Line 1 : gridx,gridy,gridz
	! Line 2 : ax,ay,az,bx,by,bz,cx,cy,cz
	! Line 3 : originx, originy, originz
	! Line 4 : loop order (e.g. 'zyx')
	! Line 5+: data (N = gridx*gridy*gridz)

	! Get file data
	write(0,*) "Reading file..."
	success = readline(11)
	if (nargsparsed.ne.3) stop "Error reading gridpoint specification."
	ngrid(1) = argi(1)
	ngrid(2) = argi(2)
	ngrid(3) = argi(3)
	success = readline(11)
	if (nargsparsed.ne.9) stop "Error reading grid axes specification."
	do n=1,9
	  axes(n) = argr(n)
	end do
	success = readline(11)
	if (nargsparsed.ne.3) stop "Error reading grid origin specification."
	origin(1) = argr(1)
	origin(2) = argr(2)
	origin(3) = argr(3)
	success = readline(11)
	if (nargsparsed.ne.1) stop "Error reading loop order specification."
	looporder = arg(1)
	do n=1,3
	  ! Grab character and convert to lowercase
	  x = ichar(looporder(n:n))
	  if (x.lt.120) x = x + 32
	  if ((x.lt.120).or.(x.gt.122)) stop "Illegal character found in loop specification."
	  loop(n) = x-119
	end do

	! Read in data
	allocate(grid(ngrid(1),ngrid(2),ngrid(3)))
	x = 1
	y = 1
	z = 1
	do n=1,ngrid(1)*ngrid(2)*ngrid(3)
	  point(loop(1)) = x
	  point(loop(2)) = y
	  point(loop(3)) = z
	  read(11,"(f20.12)") grid(point(1),point(2),point(3))
	  ! Increase counters
	  x = x + 1
	  if (x.gt.ngrid(loop(1))) then
	    x = 1
	    y = y + 1
	    if (y.gt.ngrid(loop(2))) then
	      y = 1
	      z = z + 1
	      if (z.gt.ngrid(loop(3))) write(0,*) "Array is full."
	    end if
	  end if
	end do
	close(11)

	! Print summary of input information
	write(0,"(a)") "-------------"
	write(0,"(a)") "Original Data"
	write(0,"(a)") "-------------"
	write(0,"(a,3i5)") "Gridpoints  : ",ngrid
	write(0,"(a,3f12.6)") "Grid x-axis : ",(axes(n),n=1,3)
	write(0,"(a,3f12.6)") "     y-axis : ",(axes(n),n=4,6)
	write(0,"(a,3f12.6)") "     z-axis : ",(axes(n),n=7,9)
	write(0,"(a,3f12.6)") "Grid origin : ",origin
	write(0,"(a,3i5)") "Loop order  : ",loop

	! Determine number of gridpoints to step in each direction, based on gaussian width
	extents = 1
	do n=1,3
	  mag(n) = sqrt(axes((n-1)*3+1)*axes((n-1)*3+1) + axes((n-1)*3+2)*axes((n-1)*3+2) + axes((n-1)*3+3)*axes((n-1)*3+3))
	  do while ( exp(  (-((mag(n)*extents(n))*(mag(n)*extents(n))))  /  (2*sigma*sigma)) .gt.1.0e-2 )
	    extents(n) = extents(n) + 1
	  end do
	end do
	write(0,*) "Extents for gaussian filtering (in gridpoints): ", extents

	! Create space for filtered data
	allocate(gaussgrid(ngrid(1),ngrid(2),ngrid(3)))
	gaussgrid = 0.0
	twosigma2 = 2.0 * sigma * sigma

	! Loop over original data
	do x=1,ngrid(1)
	  write(0,"(a,i4,a,i4)") "At x = ", x, " of ", ngrid(1)
	  do y=1,ngrid(2)
	    do z=1,ngrid(3)

	      ! Loop over gaussian extent
	      do i=-extents(1),extents(1)
	        do j=-extents(2),extents(2)
	          do k=-extents(3),extents(3)
		    point(1) = x+i
		    if ((point(1).lt.1).or.(point(1).gt.ngrid(1))) cycle
		    point(2) = y+j
		    if ((point(2).lt.1).or.(point(2).gt.ngrid(2))) cycle
		    point(3) = z+k
		    if ((point(3).lt.1).or.(point(3).gt.ngrid(3))) cycle

		    xyz = (mag(1)*i)*(mag(1)*i) + (mag(2)*j)*(mag(2)*j) * (mag(3)*k)*(mag(3)*k)
		    g = grid(x,y,z) * exp(-xyz/twosigma2)
		    gaussgrid(point(1),point(2),point(3)) = gaussgrid(point(1),point(2),point(3)) + g
		  end do
		end do
	      end do

	    end do
	  end do
	end do

	! Normalise
	!gaussgrid = gaussgrid / (4.0*pi*pi*sigma*sigma*sigma*sigma)
	!gaussgrid = gaussgrid * (sum(grid) / sum(gaussgrid))
	gaussgrid = gaussgrid * (maxval(grid) / maxval(gaussgrid))
	

	write(0,"(a)") ""
	write(0,"(a)") "-------------"
	write(0,"(a)") "Modified Data"
	write(0,"(a)") "-------------"
	write(0,"(a,3i5)") "Gridpoints  : ",ngrid
	write(0,"(a,3f12.6)") "Grid x-axis : ",(axes(n),n=1,3)
	write(0,"(a,3f12.6)") "     y-axis : ",(axes(n),n=4,6)
	write(0,"(a,3f12.6)") "     z-axis : ",(axes(n),n=7,9)
	write(0,"(a,3f12.6)") "Grid origin : ",origin
	write(0,"(a,3i5)") "Loop order  : ",loop

        open(unit=11,file=outfile,form='formatted',status='new')
	write(11,*) ngrid
	write(11,"(9f6.2)") axes
	write(11,"(3f10.4)") origin
	write(11,*) looporder
	
	do x=1,ngrid(1)
	  do y=1,ngrid(2)
	    do z=1,ngrid(3)
	      point(loop(3)) = x
	      point(loop(2)) = y
	      point(loop(1)) = z
	  !write(55,*) "writing point",x,y,z
	      write(11,"(f12.5)") gaussgrid(point(1),point(2),point(3))
	    end do
	  end do
	end do

	close(11)

	end program pdensgauss
