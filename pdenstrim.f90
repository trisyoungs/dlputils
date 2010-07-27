	! Trim (reduce) the data within a pdens file

	program pdenstrim
	use parse
	implicit none
	character*80 :: infile, outfile
	character*20 :: looporder, temp
	character :: c
	integer :: ngrid(3), nargs, n, x, y, z, loop(3), point(3), nchanges = 0
	integer :: mingrid(3), maxgrid(3), d(3)
	logical :: gnufile = .false., success
	real*8 :: axes(9), origin(3), v(3)
	real*8, allocatable :: dat(:,:,:)
	integer :: iargc

        nargs = iargc()
        if (nargs.lt.2) stop "Usage: pdenstrim <pdensfile> <outputfile> [-min x y z] [-max x y z] [-gnuplot]"
        call getarg(1,infile)
        call getarg(2,outfile)

	! Set defaults
	mingrid = (/ 1, 1, 1 /)
	maxgrid = (/ -1, -1, -1 /)
	n = 2
        do
          n = n + 1; if (n.gt.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-min")
              n = n + 1; call getarg(n,temp); read(temp,"(i10)") mingrid(1)
              n = n + 1; call getarg(n,temp); read(temp,"(i10)") mingrid(2)
              n = n + 1; call getarg(n,temp); read(temp,"(i10)") mingrid(3)
              write(0,"(A,3I4)") "New grid will begin at (XYZ) gridpoints : ",mingrid
	      nchanges = nchanges + 1
            case ("-max")
              n = n + 1; call getarg(n,temp); read(temp,"(i10)") maxgrid(1)
              n = n + 1; call getarg(n,temp); read(temp,"(i10)") maxgrid(2)
              n = n + 1; call getarg(n,temp); read(temp,"(i10)") maxgrid(3)
              write(0,"(A,3I4)") "New grid will begin at (XYZ) gridpoints : ",maxgrid
	      nchanges = nchanges + 1
            case ("-gnuplot")
	      gnufile = .true.
	      write(0,*) "Files for GnuPlot will be written instead."
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

	! Apply modifications (if any)
	if (nchanges.eq.0) stop

	! Set existing values if new values were not provided
	do n=1,3
	  if (mingrid(n).eq.-1) mingrid(n) = 1
	  if (maxgrid(n).eq.-1) maxgrid(n) = ngrid(n)
	end do

	! Work out the new grid origin
	d = mingrid - 1
	!write(0,*) d
	v(1) = origin(1) + d(1)*axes(1) + d(2)*axes(4) + d(3)*axes(7)
	v(2) = origin(2) + d(1)*axes(2) + d(2)*axes(5) + d(3)*axes(8)
	v(3) = origin(3) + d(1)*axes(3) + d(2)*axes(6) + d(3)*axes(9)
	origin = v

	! Write out the modified data
	x = mingrid(1)
	y = mingrid(2)
	z = mingrid(3)
	d = maxgrid - mingrid
	d = d + 1
	ngrid = d

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

	if (.not.gnufile) then

          open(unit=11,file=outfile,form='formatted',status='new')
	  write(11,*) ngrid
	  write(11,"(9f6.2)") axes
	  write(11,"(3f10.4)") origin
	  write(11,*) looporder
	
	  do x=mingrid(1),maxgrid(1)
	    do y=mingrid(2),maxgrid(2)
	      do z=mingrid(3),maxgrid(3)
	        point(loop(3)) = x
	        point(loop(2)) = y
	        point(loop(1)) = z
	  !write(55,*) "writing point",x,y,z
	        write(11,"(f12.5)") dat(point(1),point(2),point(3))
	      end do
	    end do
	  end do

	else

          open(unit=11,file=outfile,form='formatted',status='new')
	  do x=mingrid(1),maxgrid(1)
	    do y=mingrid(2),maxgrid(2)
	      do z=mingrid(3),maxgrid(3)
	        point(loop(3)) = x
	        point(loop(2)) = y
	        point(loop(1)) = z
	  !write(55,*) "writing point",x,y,z
	        write(11,"(f12.5)") dat(point(1),point(2),point(3))
	      end do
	    end do
	    write(11,*) ""
	  end do

	end if
	close(11)

	end program pdenstrim
