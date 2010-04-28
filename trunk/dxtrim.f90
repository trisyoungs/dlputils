	! Simple prog to trim (DX) grid data down to the specified size

	program dxtrim
	implicit none
	character*80 :: datafile
	character*20 :: temparg
	real*8, allocatable :: dat(:,:,:)
	integer :: iargc
	integer :: nargs, ngrid, nx, ny, nz, i1, i2, i3
	
	nargs = iargc()
	if (nargs.NE.5) stop "Usage : dxmirror <datafile> <ngrid> <newgridx> <newgridy> <newgridz>"
	call getarg(1,datafile)
	call getarg(2,temparg); read(temparg,"(i6)") ngrid
	call getarg(3,temparg); read(temparg,"(i6)") nx
	call getarg(4,temparg); read(temparg,"(i6)") ny
	call getarg(5,temparg); read(temparg,"(i6)") nz
	
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

	! Do the trim, writing out the new data
	write(0,*) "Preforming trim..."

	open(unit=11,file=datafile,form="formatted",status="replace")
	do i1 = -ngrid,ngrid
	  if (abs(i1).gt.nx) cycle
	  do i2 = -ngrid,ngrid
	    if (abs(i2).gt.ny) cycle
	    do i3 = -ngrid,ngrid
	      if (abs(i3).gt.nz) cycle
	      write(11,"(f10.5)") dat(i1,i2,i3)
	    end do
	  end do
	end do
	close(11)

	write(0,*) "Done."
	end program dxtrim
