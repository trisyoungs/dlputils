	! ###############
	! dlp filter
	! ###############

	program everydlp
	use dlprw
	implicit none
	character*80 :: hisfile,outfile,temp
	integer :: nargs,n,i,m,o,success,nframes,nfiltered,interval
	integer :: iargc

	nargs = iargc()
	if (nargs.NE.3) stop "Usage : everydlp <input HISfile> <output HISfile> <interval>"
	call getarg(1,hisfile)
	call getarg(2,outfile)
	call getarg(3,temp); read(temp,"(i5)") interval

	write(6,*) "Filtering..."
	nfiltered = 0

	nframes = 0
	! Open the file and read in the file header (not used)
	call openhis(hisfile,15)
        success = readheader()
	if (success.EQ.-1) stop "Couldn't read header of first file!"
	write(0,*) "Writing header of new file..."
	if (writeheader(outfile,14,0).EQ.-1) stop "Failed to write new history file header!"

10	success = readframe()
	if (success.NE.0) then
	  write(0,"(A,A)") "Failed to read first frame of file."
	  stop "Quit."
	end if
	nframes = nframes + 1
	if (mod(nframes,interval).eq.0) then
	  ! Write out the frame...
	  success = writeframe()
	  if (success.EQ.-1) stop "Failed to write frame!"
	  nfiltered = nfiltered + 1
	end if
	goto 10

14	write(0,*) "Found end of frames file."
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Filtered ",nfiltered," frames from file."
	goto 20

17	stop "Error reading from frames file."
	  
20      close(14)
	stop "Finished filtering."

	end program everydlp
