	! ###############
	! dlp reorder
	! ###############

	program dlpreorder
	use dlprw
	implicit none
	character*80 :: hisfile,temp,outfile,newhisfile
	real*8, allocatable :: newx(:), newy(:), newz(:), newfx(:), newfy(:), newfz(:), newvx(:), newvy(:), newvz(:)
	real*8, allocatable :: newmass(:), newcharge(:)
	character*8, allocatable :: newatmname(:)
	integer :: nargs,n,i,m,o,success,nframes, sp, neworder(100), newnatoms, oldnatoms
	integer :: iargc
	logical :: failed_header

	nargs = iargc()
	if (nargs.le.2) stop " Usage: dlpreorder <HISTORYfile> <OUTPUTfile> <sp1> [sp2] [sp3]..."
	call getarg(1,hisfile)
	call getarg(2,outfile)
	neworder = 0
	do n=3,nargs
	  call getarg(n,temp)
	  read(temp, "(i3)") neworder(n-2)
	end do

	newhisfile="reordered.HISu"

	nframes = 0
	if (outinfo(outfile,1).eq.-1) stop "Couldn't read OUT file."
	! Establish new number of atoms and check species limit
	newnatoms = 0
	do n=1,100
	  if (neworder(n).gt.nspecies) stop "Invalid species number given in list."
	  newnatoms = newnatoms + s_nmols(neworder(n)) * s_natoms(neworder(n))
	end do

	! Open the file and read in the file header
	call openhis(hisfile,15)
        success = readheader()
	if (success.EQ.-1) stop "Couldn't read header of first file!"
	if (newnatoms.gt.natms) stop "New number of atoms cannot be greater than the old number"
	write(0,*) "Writing header of new file..."
	oldnatoms = natms
	natms = newnatoms
	if (writeheader(newhisfile,14,0).EQ.-1) stop "Failed to write new history file header!"
	natms = oldnatoms

	! Allocate arrays
	allocate (newx(natms))
	allocate (newy(natms))
	allocate (newz(natms))
	allocate (newfx(natms))
	allocate (newfy(natms))
	allocate (newfz(natms))
	allocate (newvx(natms))
	allocate (newvy(natms))
	allocate (newvz(natms))
	allocate (newmass(natms))
	allocate (newcharge(natms))
	allocate (newatmname(natms))

	write(6,*) "Reordering..."

10	success = readframe()
	if (success.NE.0) goto 15
	nframes = nframes + 1
	if (MOD(nframes,100).EQ.0) write(0,*) nframes
	! Reorder atom data
	newx = 0.0
	newy = 0.0
	newz = 0.0
	newfx = 0.0
	newfy = 0.0
	newfz = 0.0
	newvx = 0.0
	newvy = 0.0
	newvz = 0.0
	i = 0
	do n=1,100
	  sp = neworder(n)
	  if (sp.eq.0) cycle
	  do m=1,s_nmols(sp)
	    do o=0,s_natoms(sp)-1
	      i = i + 1
	      newx(i) = xpos(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newy(i) = ypos(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newz(i) = zpos(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newfx(i) = xfor(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newfy(i) = yfor(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newfz(i) = zfor(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newvx(i) = xvel(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newvy(i) = yvel(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newvz(i) = zvel(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newatmname(i) = atmname(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newcharge(i) = charge(s_start(sp)+o+(m-1)*s_natoms(sp))
	      newmass(i) = mass(s_start(sp)+o+(m-1)*s_natoms(sp))
	    end do
	  end do
	end do

	! Store new data in old arrays
	do i=1,newnatoms
	  xpos(i) = newx(i)
	  ypos(i) = newy(i)
	  zpos(i) = newz(i)
	  xfor(i) = newfx(i)
	  yfor(i) = newfy(i)
	  zfor(i) = newfz(i)
	  xvel(i) = newvx(i)
	  yvel(i) = newvy(i)
	  zvel(i) = newvz(i)
	end do

	! Write out the frame...
	natms = newnatoms
	success = writeframe()
	natms = oldnatoms
	if (success.EQ.-1) stop "Failed to write frame!"
	goto 10
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Read/wrote ",nframes," from file ",hisfile
	  
20      close(14)
	stop "Finished concatenation."

	end program dlpreorder
