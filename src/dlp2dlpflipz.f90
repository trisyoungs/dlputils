	! ###############
	! dlp converter - and flip z coordinate about zero
	! ###############

	program dlp2flipx
	use dlprw
	implicit none
	character*80 :: hisfile,temp,outfile
	integer :: nargs,n,i,m,success,nframes,destfmt
	integer :: iargc
	logical :: failed_header

	nargs = iargc()
	if (nargs.NE.2) stop " Usage: dlp2 <HISTORYfile> <format (0=unformatted, 1=formatted)"
	call getarg(1,hisfile)
	call getarg(2,temp)
	read(temp,"(I2)") destfmt

	write(6,*) "Converting..."
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
	! Swap z coordinate sign (and velocity and force too)
	do n=1,natms
	  zpos(n) = -zpos(n)
	  zvel(n) = -zvel(n)
	  zfor(n) = -zfor(n)
	end do
	! Write out the frame...
	success = writeframe()
	if (success.EQ.-1) stop "Failed to write frame!"
	goto 10
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Read/wrote ",nframes," from file ",hisfile
	  
20      close(14)
	stop "Finished concatenation."

	end program dlp2flipx
