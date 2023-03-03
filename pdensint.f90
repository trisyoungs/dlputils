	! Integrate pdens file for various cutoffs

	program pdensint
	use parse; use pdensrw ; use utility
	implicit none
	character*80 :: infile, outfile
	character*20 :: looporder, temp
	character :: c
	integer :: nargs, n, m, x, y, z, point(3), nsteps = 10000
	real*8 :: v(3), maximum, nmols, delta, cutoff, vol, total, tgt, tgtdelta
	real*8, allocatable :: ints(:,:)
	type(PDens) :: grid
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
	if (.not.loadPDens(infile,grid)) stop "Couldn't load specified pdens file."

	maximum = maxval(grid%grid)
	total = sum(grid%grid)

	! Calculate volume integrations for various cutoffs
	vol = volume(grid%axes)
	delta = maximum / nsteps
	total = total * vol
	write(0,"(a,e12.5)") "Maximum value found in data is : ", maximum
	write(0,"(a,f10.6)") "Cubic volume of grid voxel : ", vol
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
	  do x=grid%gridMin(1),grid%gridMax(1)
	    do y=grid%gridMin(2),grid%gridMax(2)
	      do z=grid%gridMin(3),grid%gridMax(3)
	        if (grid%grid(x,y,z) > cutoff) nmols = nmols + grid%grid(x,y,z)
	      end do
	    end do
	  end do
	  nmols = nmols * vol
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
