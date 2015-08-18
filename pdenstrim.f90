	! Trim (reduce) the data within a pdens file

	program pdenstrim
	use parse; use PDensRW
	implicit none
	character*80 :: infile, outfile
	character*20 :: temp
	character :: c
	type(PDens) :: original, trimmed
	integer :: ngrid(3), nargs, n, i, j, k, point(3), nchanges = 0
	integer :: mingrid(3), maxgrid(3), d(3)
	logical :: gnufile = .false., success, minProvided = .false., maxProvided = .false.
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
              n = n + 1; mingrid(1) = getargi(n)
              n = n + 1; mingrid(2) = getargi(n)
              n = n + 1; mingrid(3) = getargi(n)
              write(0,"(A,3I4)") "New grid will begin at (XYZ) gridpoints : ",mingrid
	      nchanges = nchanges + 1
	      minProvided = .true.
            case ("-max")
              n = n + 1; maxgrid(1) = getargi(n)
              n = n + 1; maxgrid(2) = getargi(n)
              n = n + 1; maxgrid(3) = getargi(n)
              write(0,"(A,3I4)") "New grid will end at (XYZ) gridpoints : ",maxgrid
	      nchanges = nchanges + 1
	      maxProvided = .true.
            case ("-gnuplot")
	      gnufile = .true.
	      write(0,*) "Files for GnuPlot will be written instead."
	    case default
	      write(0,*) "Unrecognised argument :",temp
	      stop
	  end select
	end do

	! Open existing pdens
	if (.not.loadPDens(infile, original)) stop "Failed to load original pdens."

	if (nchanges.eq.0) stop

	! Set existing values if new values were not provided
	if (.not.minProvided) mingrid = original%gridMin
	if (.not.maxProvided) maxgrid = original%gridMax

	! Work out the new grid origin
	d = mingrid - 1
	trimmed%origin(1) = original%origin(1) + d(1)*original%axes(1) + d(2)*original%axes(4) + d(3)*original%axes(7)
	trimmed%origin(2) = original%origin(2) + d(1)*original%axes(2) + d(2)*original%axes(5) + d(3)*original%axes(8)
	trimmed%origin(3) = original%origin(3) + d(1)*original%axes(3) + d(2)*original%axes(6) + d(3)*original%axes(9)

	! Copy axes and loop from old pdens
	trimmed%axes = original%axes
	trimmed%loop = original%loop

	! Determine new grid extents
	d = maxgrid - mingrid
	d = d + 1
	ngrid = d

	! Allocate new pdens space
	if (.not.allocPDens(trimmed,mingrid(1),mingrid(2),mingrid(3),maxgrid(1),maxgrid(2),maxgrid(3))) stop
	do i=mingrid(1),maxgrid(1)
	  do j=mingrid(2),maxgrid(2)
	    do k=mingrid(3),maxgrid(3)
	      trimmed%grid(i,j,k) = original%grid(i,j,k)
	    end do
	  end do
	end do

	if (.not.gnufile) then

	  if (.not.savePDens(outfile, trimmed)) stop

	else

          open(unit=11,file=outfile,form='formatted',status='new')
	  do i=mingrid(original%loop(1)),maxgrid(original%loop(1))
	    do j=mingrid(original%loop(2)),maxgrid(original%loop(2))
	      do k=mingrid(original%loop(3)),maxgrid(original%loop(3))
	        point(original%loop(1)) = i
	        point(original%loop(2)) = j
	        point(original%loop(3)) = k
	  !write(55,*) "writing point",x,y,z
	        write(11,"(f12.5)") trimmed%grid(point(1),point(2),point(3))
	      end do
	      write(11,*) ""
	    end do
	  end do
	  close(11)

	end if

	end program pdenstrim
