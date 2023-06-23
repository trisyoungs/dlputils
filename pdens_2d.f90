	! Analysis of 3d pdens file, to obtain 2d slices and 2d integration
	! Updated 24/08/2018 to actually write to a file

	program pdens_2d
	use parse; use pdensrw ; use utility
	implicit none
	character*80 :: infile, outfile, basename, resfile
	character*20 :: looporder, temp
	character(len=2) :: slice
	integer :: nargs, n, m, x, y, z, point(3), mid, baselen
	real*8 :: v(3), maximum, nmols, delta, cutoff, vol, total, tgt, tgtdelta
	real*8, allocatable :: ints(:,:)
	type(PDens) :: grid
	integer :: iargc
	logical :: success, inflect
	logical :: intr=.false.

        nargs = iargc()
        if (nargs.lt.1) stop "Usage: pdens_2d <pdensfile> [-slice xy|xz|yz]"
        call getarg(1,infile)

	! Get CLI arguments
	n = 1
        do
	n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-slice")
	      n = n + 1; call getarg(n,temp); read(temp,"(A2)") slice
	    if (slice.EQ."xy".OR.slice.EQ."xz".OR.slice.EQ."yz") then
	    	write(0,*) "slice chosen = ", slice
		else
		write(0,*) temp, "unrecognised option for 2d slice, choose either xy|xz|yz"
		stop
	    endif
	    !case ("-int")
		!intr = .true.
	      !write(0,*) "2d slices based on integral over all of other axes [NOT CURRENTLY WORKING]"
		!stop		
	    case default
	      write(0,*) "Unrecognised argument :",temp
	      stop
	  end select
	end do

        ! Input file
	write(0,"(A,A)") "Input .pdens file : ",infile
	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  IF (infile(n:n).eq.".") then
	    baselen=n
	    goto 50
	  endIF
	end do
50	if (baselen.eq.-1) then
	  basename="2ddistresults."
	  baselen=14
	else
	  basename=infile(1:baselen)
	endif

	!open output file for writing
	resfile=basename(1:baselen)//slice//".pdens2d"
	open(unit=1,file=resfile,form="formatted")


	! Open file
	if (.not.loadPDens(infile,grid)) stop "Couldn't load specified pdens file."

	! Print xy plane e.g. z=0
	! mid point of xy
	!mid = grid%gridMax(3) + 1

	if (slice.EQ."xy")then

	do m=grid%gridMin(1),grid%gridMax(1)
	   do n=grid%gridMin(2),grid%gridMax(2)
		write (1,"(3(F12.8,1x))") (m*grid%axes(1)), (n*grid%axes(5)), grid%grid(m,n,0)
	   enddo
	write(1,*) " "
        enddo

	elseif (slice.EQ."xz")then

	do m=grid%gridMin(1),grid%gridMax(1)
	   do n=grid%gridMin(3),grid%gridMax(3)
		write (1,"(3(F12.8,1x))") (m*grid%axes(1)), (n*grid%axes(9)), grid%grid(m,0,n)
	   enddo
	write(1,*) " "
        enddo

	elseif (slice.EQ."yz") then

	do m=grid%gridMin(2),grid%gridMax(2)
	   do n=grid%gridMin(3),grid%gridMax(3)
		write (1,"(3(F12.8,1x))") (m*grid%axes(5)), (n*grid%axes(9)), grid%grid(0,m,n)
	   enddo
	write(1,*) " "
        enddo

	!else write(*,*) "error unrecognised option for -slice"
	endif

	end program pdens_2d
