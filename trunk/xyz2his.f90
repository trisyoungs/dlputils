	! #####################
	! XYZ to HISu converter
	! #####################

	program xyz2his
	use dlprw; use parse
	implicit none
	character*80 :: xyzfile,temp,outfile,header
	integer :: nargs,n,i,m,success,nframes,destfmt
	integer :: iargc
	logical :: hascell = .FALSE.

	nargs = iargc()
	if (nargs.lt.2) stop " Usage: xyz2his <xyz file> <output.HISu> [ax ay az bx by bz cx cy cz]"
	call getarg(1,xyzfile)
	call getarg(2,outfile)
	imcon = 0
	if (nargs.gt.2) then
	  if (nargs.eq.11) then
	    do n=1,9
	      call getarg(n+2,temp)
	      read(temp, "(f20.14)") cell(n)
	    end do
	    hascell = .TRUE.
	    imcon = 2
	  else
	    stop "Error: Partial cell specification given?"
	  end if
	end if

	! Set destination format (0 = unformatted history file, 1 = formatted history file)
	destfmt = 0
	
	! Open the xyz file and read in the first frame to get atom numbers and names
	open(unit=15,file=xyzfile,form='formatted',status='old',err=999)
	read(15,"(i20)") natms
	call alloc_xyz
	read(15,"(a80)") header
	do n=1,natms
	  if (.not.readline(15)) goto 999
	  atmname(n) = arg(1)
	  keytrj = 0
	  if (nargsparsed.ge.7) keytrj = 1
	  if (nargsparsed.eq.10) keytrj = 2
	end do
	rewind(15)

	nframes = 0
	write(0,*) "Writing header of new file..."
	if (writeheader(outfile,14,destfmt).EQ.-1) stop "Failed to write new history file header!"

10	read(15,"(i20)",end=15,err=15) natms
	read(15,"(a80)") header
	do n=1,natms
	  if (.not.readline(15)) goto 999
	  xpos(n) = argr(2)
	  ypos(n) = argr(3)
	  zpos(n) = argr(4)
	  if (nargsparsed.ge.7) then
	    xvel(n) = argr(5)
	    yvel(n) = argr(6)
	    zvel(n) = argr(7)
	  end if
	  if (nargsparsed.eq.10) then
	    xfor(n) = argr(8)
	    yfor(n) = argr(9)
	    zfor(n) = argr(10)
	  end if
	end do
	nframes = nframes + 1
	if (MOD(nframes,100).EQ.0) write(0,*) nframes
	! Write out the frame...
	success = writeframe()
	if (success.EQ.-1) stop "Failed to write frame!"
	goto 10
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Read/wrote ",nframes," from file ",xyzfile
	  
20      close(14)
	stop "Finished conversion."
999	stop "Something bad happened!"

	end program xyz2his
