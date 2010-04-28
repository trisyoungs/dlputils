	! Integrate pdens file for various cutoffs

	program pdensint
	use parse
	implicit none
	character*80 :: infile, outfile
	character*20 :: looporder, temp
	character :: c
	integer :: ngrid(3), nargs, n, m, x, y, z, loop(3), point(3), nsteps = 10000
	real*8 :: axes(9), origin(3), v(3), maximum, nmols, delta, cutoff, volume, total, tgt, tgtdelta
	real*8, allocatable :: dat(:,:,:), ints(:,:)
	integer :: iargc
	logical :: success, inflect

        nargs = iargc()
        if (nargs.lt.1) stop "Usage: pdensperc <pdensfile> [-nsteps n]"
        call getarg(1,infile)

	! Get CLI arguments
	n = 1
        do
          n = n + 1; if (n.gt.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-nsteps")
	      n = n + 1; call getarg(n,temp); read(temp,"(I4)") nsteps
	    case default
	      write(0,*) "Unrecognised argument :",temp
	      stop
	  end select
	end do

	write(0,*) "There will be ", nsteps, "points calculated."
	allocate(ints(nsteps,3))
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
	maximum = -10000.0
	total = 0.0
	allocate(dat(ngrid(1),ngrid(2),ngrid(3)))
	x = 1
	y = 1
	z = 1
	do n=1,ngrid(1)*ngrid(2)*ngrid(3)
	  point(loop(1)) = x
	  point(loop(2)) = y
	  point(loop(3)) = z
	  read(11,"(f20.12)") dat(point(1),point(2),point(3))
	  total = total + dat(point(1),point(2),point(3))
	  if (dat(point(1),point(2),point(3)).gt.maximum) maximum = dat(point(1),point(2),point(3))
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

	! Calculate volume integrations for various cutoffs
	volume = axes(1)*axes(5)*axes(9)
	delta = maximum / nsteps
	total = total * volume
	write(0,"(a,e12.5)") "Maximum value found in data is : ", maximum
	write(0,"(a,f10.6)") "Cubic volume of grid voxel : ", volume
	write(0,"(a,f10.6)") "Total integral of grid : ", total
	write(0,"(a,f10.6)") "Delta is : ", delta
	cutoff = 0.0

	n = 0
	write(6,"(a)") "# Cutoff        NMols           %Total"
	do
	  n = n + 1
	  if (cutoff.gt.maximum) exit
	  ! Calculate
	  nmols = 0.0
	  do x=1,ngrid(1)
	    do y=1,ngrid(2)
	      do z=1,ngrid(3)
	        if (dat(x,y,z) > cutoff) nmols = nmols + dat(x,y,z)
	      end do
	    end do
	  end do
	  nmols = nmols * volume
	  ints(n,1:3) = (/ cutoff, nmols, nmols / total * 100.0 /)
	  write(6,"(f12.6,4x,f12.6,4x,f12.6)") (ints(n,m),m=1,3)
	  cutoff = cutoff + delta
	end do

	! Find significant values
	tgt = 5.0
	tgtdelta = 100.0
	do n=nsteps,1,-1
	  !write(0,*) n,ints(n,3)-tgt
	  if (dabs(ints(n,3)-tgt).le.tgtdelta) then
	    tgtdelta = dabs(ints(n,3)-tgt)
	  else
	    write(0,"(a,f4.1,a,f8.5,a,f10.5,a)") "For top ",tgt,"% of molecules a cutoff of ", ints(n,1), " is required (", ints(n,2)," molecules)"
	    tgt = tgt * 2.0
	    if (tgt.gt.100) exit
	    tgtdelta = 100.0
	  end if
	end do

	end program pdensint
