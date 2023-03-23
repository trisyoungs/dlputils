	! Simple prog to symmetrize pdens data in the given direction (about the centrepoint of the data)

	program pdensmirror
	use parse
	implicit none
	character*80 :: datafile
	character*20 :: temparg
	character*3 :: direction, looporder
	real*8, allocatable :: dat(:,:,:)
	integer :: nargs, ngrid(3), midgrid(3), n, i1, i2, i3, x, y, z, c, point(3), loop(3)
	integer :: iargc
	real*8 :: origin(3), axes(9), avg
	logical :: success
	
	nargs = iargc()
	if (nargs.ne.2) stop "Usage : pdensmirror <datafile> <direction(s)>"
	call getarg(1,datafile)
	call getarg(2,direction)
	
	! Open file
	open(unit=11, file=datafile, form='formatted', status='old')

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
	allocate(dat(ngrid(1),ngrid(2),ngrid(3)))
	x = 1
	y = 1
	z = 1
	do n=1,ngrid(1)*ngrid(2)*ngrid(3)
	  point(loop(1)) = x
	  point(loop(2)) = y
	  point(loop(3)) = z
	  read(11,"(f20.12)") dat(point(1),point(2),point(3))
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

	midgrid = (ngrid-1)/2 + 1
	write(0,*) midgrid

	! Do the mirroring
	do c=1,3
	  if (direction(c:c).eq." ") exit
	  write(0,"('Mirroring in ',a1,' direction...')") direction(c:c)

	  select case (direction(c:c))
	    case ("x")
	      do i1 = 1,midgrid(1)
		do i2 = 1, ngrid(2)
		  do i3 = 1, ngrid(3)
		    ! Take the average of the current value with it's mirror and put it back into the array
		    avg = (dat(i1,i2,i3) + dat(ngrid(1)-(i1-1),i2,i3)) * 0.5
		    dat(i1,i2,i3) = avg
		    dat(ngrid(1)-(i1-1),i2,i3) = avg
		  end do
		end do
	      end do
	    case ("y")
	      do i1 = 1,ngrid(1)
		do i2 = 1, midgrid(2)
		  do i3 = 1, ngrid(3)
		    ! Take the average of the current value with it's mirror and put it back into the array
		    avg = (dat(i1,i2,i3) + dat(i1,ngrid(2)-(i2-1),i3)) * 0.5
		    dat(i1,i2,i3) = avg
		    dat(i1,ngrid(2)-(i2-1),i3) = avg
		  end do
		end do
	      end do
	    case ("z")
	      do i1 = 1,ngrid(1)
		do i2 = 1, ngrid(2)
		  do i3 = 1, midgrid(3)
		    ! Take the average of the current value with it's mirror and put it back into the array
		    avg = (dat(i1,i2,i3) + dat(i1,i2,ngrid(3)-(i3-1))) * 0.5
		    dat(i1,i2,i3) = avg
		    dat(i1,i2,ngrid(3)-(i3-1)) = avg
		  end do
		end do
	      end do
	  end select
	end do

	! Write new data
        open(unit=11,file=datafile,form='formatted',status='replace')
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
	      write(11,"(f12.5)") dat(point(1),point(2),point(3))
	    end do
	  end do
	end do
	close(11)


	write(0,*) "Done."
	end program pdensmirror
