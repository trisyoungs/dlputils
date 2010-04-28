	! Simple prog to symmetrize (DX) grid data in the given direction (about the centre of the data)

	program dxmirror
	implicit none
	character*80 :: datafile
	character*20 :: temparg
	character*3 :: direction
	real*8, allocatable :: dat(:,:,:)
	integer :: nargs, ngrid, n1, n2, n3, m1, m2, m3, i1, i2, i3, c
	integer :: iargc
	real*8 :: mag, avg
	
	nargs = iargc()
	if (nargs.NE.3) stop "Usage : dxmirror <datafile> <ngrid> <direction>"
	call getarg(1,datafile)
	call getarg(2,temparg); read(temparg,"(i6)") ngrid
	call getarg(3,direction)
	
	allocate(dat(-ngrid:ngrid,-ngrid:ngrid,-ngrid:ngrid))

	! Read in the data
	write(0,*) "Reading original data..."
	open(unit=10,file=datafile,form="formatted",status="old")
	do i1 = -ngrid,ngrid
	  do i2 = -ngrid,ngrid
	    do i3 = -ngrid,ngrid
	      read(10,"(f20.12)") dat(i1,i2,i3)
	    end do
	  end do 
	end do
	close(10)

	! Do the mirroring
	do c=1,3
	if (direction(c:c).eq." ") exit
	write(0,"('Mirroring in ',a1,' direction...')") direction(c:c)

	do i1 = -ngrid,ngrid
	  do i2 = -ngrid,ngrid
	    do i3 = 0,ngrid
	      ! The first two loops correspond to the directions we explore fully.
	      ! The last loop corresponds to the mirror direction
	      select case (direction(c:c))
		case ("x")
		  n1 = i3; m1 = -i3
		  n2 = i1; m2 = i1
		  n3 = i2; m3 = i2
		case ("y")
		  n1 = i1; m1 = i1
		  n2 = i3; m2 = -i3
		  n3 = i2; m3 = i2
		case ("z")
		  n1 = i1; m1 = i1
		  n2 = i2; m2 = i2
		  n3 = i3; m3 = -i3
	      end select
	      ! Take the average of the current value with it's mirror and put it back into the array
	      avg = (dat(n1,n2,n3) + dat(m1,m2,m3)) * 0.5
	      dat(n1,n2,n3) = avg
	      dat(m1,m2,m3) = avg
	    end do
	  end do
	end do

	end do

	! Write out new data
	write(0,*) "Writing new data file..."
	open(unit=11,file=datafile,form="formatted",status="replace")
	do i1 = -ngrid,ngrid
	  do i2 = -ngrid,ngrid
	    do i3 = -ngrid,ngrid
	      write(11,"(f10.5)") dat(i1,i2,i3)
	    end do
	  end do 
	end do
	close(11)

	write(0,*) "Done."
	end program dxmirror
