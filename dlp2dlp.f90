	! ###############
	! dlp converter
	! ###############

	program dlp2dlp
	use dlprw
	implicit none
	character*80 :: hisfile,temp,outfile
	integer :: nargs,n,i,m,success,nframes,destfmt,framestodo = -1
	integer :: iargc
	logical :: failed_header

	nargs = iargc()
	if (nargs.lt.2) stop " Usage: dlp2dlp <HISTORYfile> <format (0=unformatted, 1=formatted)> [nframes]"
	call getarg(1,hisfile)
	call getarg(2,temp)
	read(temp,"(I2)") destfmt
	if (nargs.eq.3) then
	  call getarg(3,temp)
	  read(temp,"(I10)")framestodo 
	end if

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
	! Write out the frame...
	success = writeframe()
	if (success.EQ.-1) stop "Failed to write frame!"
	if (nframes.eq.framestodo) goto 15
	goto 10
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Read/wrote ",nframes," from file ",hisfile
	  
20      close(14)
	stop "Finished conversion."

	end program dlp2dlp
