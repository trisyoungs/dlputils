	! ###############
	! dlp history file converter - shift along single cartesian cell axis
	! ###############

	program dlp2dlpshift
	use dlprw
	implicit none
	character*80 :: hisfile,temp,outfile
	character*1 :: axis
	integer :: nargs,n,i,m,success,nframes,destfmt
	integer :: iargc 
	real*8 :: delta
	logical :: failed_header

	nargs = iargc()
	if (nargs.ne.4) stop " Usage: dlp2dlpshift <HISTORYfile> <format (0=unformatted, 1=formatted)> <axis (xyz)> <delta 0-1.0>"
	call getarg(1,hisfile)
	call getarg(2,temp)
	read(temp,"(I2)") destfmt
	call getarg(3,axis)
	call getarg(4,temp)
	read(temp,"(f10.6)") delta

	write(6,"(a,f9.5,a,a)") "Converting with shift ",delta," along ", axis
	if (destfmt.EQ.0) then
	  outfile="converted.HISu"
	else
	  outfile="converted.HISf"
	end if

	nframes = 0
	! Open the file and read in the file header
	call openhis(hisfile,15)
        success = readheader()
	if (success.EQ.-1) stop "Couldn't read header of first file!"
	write(0,*) "Writing header of new file..."
	if (writeheader(outfile,14,destfmt).EQ.-1) stop "Failed to write new history file header!"

10	success = readframe()
	if (success.NE.0) goto 15
	nframes = nframes + 1
	if (MOD(nframes,100).EQ.0) write(0,*) nframes

	! Shift coordinates
	if (axis(1:1).eq."x") then
	  do n=1,natms
	    xpos(n) = xpos(n) + delta * cell(1)
	    if (xpos(n).gt.cell(1)*0.5) xpos(n) = xpos(n) - cell(1)
	    if (xpos(n).lt.-cell(1)*0.5) xpos(n) = xpos(n) + cell(1)
	  end do
	else if (axis(1:1).eq."y") then
	  do n=1,natms
	    ypos(n) = ypos(n) + delta * cell(5)
	    if (ypos(n).gt.cell(5)*0.5) ypos(n) = ypos(n) - cell(5)
	    if (ypos(n).lt.-cell(5)*0.5) ypos(n) = ypos(n) + cell(5)
	  end do
	else if (axis(1:1).eq."z") then
	  do n=1,natms
	    zpos(n) = zpos(n) + delta * cell(9)
	    if (zpos(n).gt.cell(9)*0.5) zpos(n) = zpos(n) - cell(9)
	    if (zpos(n).lt.-cell(9)*0.5) zpos(n) = zpos(n) + cell(9)
	  end do
	else
	  stop "Unrecognised character given for axis."
	end if

	! Write out the frame...
	success = writeframe()
	if (success.EQ.-1) stop "Failed to write frame!"
	goto 10
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Read/wrote ",nframes," from file ",hisfile
	  
20      close(14)
	stop "Finished concatenation."

	end program dlp2dlpshift
