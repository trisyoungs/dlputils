	! Module to read/write pdens grid files

	module PDensRW

	  ! Pdens Grid Structure
	  type PDens
	    ! Number of gridpoints in each +/- direction
	    integer :: gridMin(3) = 0, gridMax(3) = 0
	    ! Data order - fastest varying (innermost) in loop first
	    integer :: loop(3) = (/ 3, 2, 1/)
	    ! Axes definition (for single voxel)
	    real*8 :: axes(9) = (/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /)
	    ! Coordinates at origin (llc) of grid
	    real*8 :: origin(3) = 0.0
	    ! Voxel grid
	    real*8, allocatable :: grid(:,:,:) 
	  end type PDens

	contains
	  
	! Clear specified structure
	subroutine clearPDens(p)
	implicit none
	type(PDens), intent(inout) :: p
	p%gridMin = 0
	p%gridMax = 0
	p%loop = (/ 3, 2, 1 /)
	p%axes(1:9) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)
	p%origin = 0.0
	if (allocated(p%grid)) deallocate(p%grid)
	end subroutine clearPDens

	! Initialise voxel array in specified structure
	logical function allocPDens(p, xmin, ymin, zmin, xmax, ymax, zmax)
	implicit none
	type(PDens), intent(inout) :: p
	integer, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
	integer :: i
	allocPDens = .true.
	if (allocated(p%grid)) deallocate(p%grid)
	p%gridMin = (/ xmin, ymin, zmin /)
	p%gridMax = (/ xmax, ymax, zmax /)
	write(0,"(a,6i4)") "Allocating PDens with min/max ", p%gridMin, p%gridMax
	allocate(p%grid(p%gridMin(1):p%gridMax(1), p%gridMin(2):p%gridMax(2), p%gridMin(3):p%gridMax(3)), stat=i)
	if (i.gt.0) then
	  write(0,*) "Failed to allocate voxel array."
	  allocPDens = .false.
	end if
	p%grid = 0.0
	end function allocPDens

	! Load specified pdens file into structure provided
	logical function loadPDens(filename, p)
	use parse
	implicit none
	character*80, intent(in) :: filename
	type(PDens), intent(out) :: p
	character*80 :: looporder
	integer :: point(3), n, i, j, k, nPoints
	logical :: success

	! Clear structure
	call clearPDens(p)

	! Attempt to open file
	open(unit=11,file=filename, form='formatted', status='old', err=999)

	! Format of data is :
	! Line 1 : ngridx,ngridy,ngridz, minx,miny,minz, maxx, maxy, maxz
	! Line 2 : ax,ay,az,bx,by,bz,cx,cy,cz
	! Line 3 : originx, originy, originz
	! Line 4 : loop order (e.g. 'zyx')
	! Line 5+: data (N = gridx*gridy*gridz)
	! Note:  Line 1 contains 3 arguments in old version, or 9 arguments in new (first three are identical)

	! Get grid extents
	write(0,*) "Reading PDens file: ", filename
	success = readline(11)
	if (nargsparsed.eq.3) then
	  write(0,*) "Old-style gridpoint specification found (2*n+1)"
	  n = (argi(1)-1)/2
	  p%gridMin = (/ -n, -n, -n /) 
	  p%gridMax = (/ n, n, n /)
	else if (nargsparsed.eq.9) then
	  write(0,*) "New-style gridpoint specification found (min,max)"
	  p%gridMin = (/ argi(4), argi(5), argi(6) /)
	  p%gridMax = (/ argi(7), argi(8), argi(9) /)
	else
	  write(0,*) "Error reading gridpoint specification."
	  loadPDens = .false.
	  return
	end if

	! Get voxel axes
	success = readline(11)
	if (nargsparsed.ne.9) then
	  write(0,*) "Error reading grid axes specification."
	  loadPDens = .false.
	  return
	end if
	do n=1,9
	  p%axes(n) = argr(n)
	end do

	! Get grid origin
	success = readline(11)
	if (nargsparsed.ne.3) then
	  write(0,*) "Error reading grid origin specification."
	  loadPDens = .false.
	  return
	end if
	p%origin(1) = argr(1)
	p%origin(2) = argr(2)
	p%origin(3) = argr(3)

	! Get loop order
	success = readline(11)
	if (nargsparsed.ne.1) then
	  write(0,*) "Error reading loop order specification."
	  loadPDens = .false.
	  return
	end if
	looporder = arg(1)
	do n=1,3
	  ! Grab character and convert to lowercase
	  i = ichar(looporder(n:n))
	  if (i.lt.120) i = i + 32
	  if ((i.lt.120).or.(i.gt.122)) then
	    write(0,*) "Illegal character found in loop specification."
	    loadPDens = .false.
	    return
	  end if
	  p%loop(n) = i-119
	end do

	! Allocate voxel array
	if (.not.allocPDens(p,p%gridMin(1),p%gridMin(2),p%gridMin(3),p%gridMax(1),p%gridMax(2),p%gridMax(3))) then
	  write(0,*) "Failed to completely load pdens file."
	  loadPDens = .false.
	  return
	end if

	! Read in voxel data
	do i=p%gridMin(p%loop(3)),p%gridMax(p%loop(3))
	  do j=p%gridMin(p%loop(2)),p%gridMax(p%loop(2))
	    do k=p%gridMin(p%loop(1)),p%gridMax(p%loop(1))
	      point(p%loop(3)) = i
	      point(p%loop(2)) = j
	      point(p%loop(1)) = k
	      read(11,*) p%grid(point(1),point(2),point(3))
	    end do
	  end do
	end do
	close(11)

	! Print summary of input information
	write(0,"(a)") "-------------"
	write(0,"(a,a)") "File: ", filename
	write(0,"(a)") "-------------"
	write(0,"(a,6i5)") "Gridmin/max  : ", p%gridMin, p%gridMax
	write(0,"(a,3f12.6)") "Grid x-axis : ", (p%axes(n),n=1,3)
	write(0,"(a,3f12.6)") "     y-axis : ", (p%axes(n),n=4,6)
	write(0,"(a,3f12.6)") "     z-axis : ", (p%axes(n),n=7,9)
	write(0,"(a,3f12.6)") "Grid origin : ", p%origin
	write(0,"(a,3i5)") "Loop order  : ", p%loop

	loadPDens = .true.
	return
998	loadPDens = .false.
	write(0,*) "Not enough data in pdens file."
	return
999	loadPDens = .false.
	write(0,*) "Failed to open file: ", filename
	return
	end function loadPDens

	! Save specified pdens structure into file
	logical function savePDens(filename, p)
	use parse
	implicit none
	character*80, intent(in) :: filename
	type(PDens), intent(in) :: p
	integer :: point(3), n, i, j, k
	logical :: success

	! Open file for writing
	open(unit=12,file=filename, form='formatted', status='replace', err=999)

	! Write description of grid / voxels
	write(12,"(9i5)") (p%gridMax(1)-p%gridMin(1))+1, (p%gridMax(2)-p%gridMin(2))+1, (p%gridMax(3)-p%gridMin(3))+1, p%gridMin, p%gridMax
	write(12,"(9(e9.3,1x))") p%axes
	write(12,"(3f10.4)") p%origin
	write(12,"(3a1)") (char(p%loop(n)+119),n=1,3)
	
	do i=p%gridMin(p%loop(3)),p%gridMax(p%loop(3))
	  do j=p%gridMin(p%loop(2)),p%gridMax(p%loop(2))
	    do k=p%gridMin(p%loop(1)),p%gridMax(p%loop(1))
	      point(p%loop(3)) = i
	      point(p%loop(2)) = j
	      point(p%loop(1)) = k
	      write(12,*) p%grid(point(1),point(2),point(3))
	    end do
	  end do
	end do
	close(12)

	savePDens = .true.
	return
999	savePDens = .false.
	write(0,*) "Failed to open file for writing: ", filename
	return
	end function savePDens

	end module PDensRW
