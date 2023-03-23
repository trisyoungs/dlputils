	! ###############
	! Add header to a restarted dlp trajectory file
	! ###############

	program addheader
	use dlprw
	implicit none
	character*80 :: headerfile,framesfile,outfile
	integer :: nargs,n,i,m,o,success,nframes
	integer :: iargc
	logical :: failed_header

	nargs = iargc()
	if (nargs.NE.3) stop " Usage <sourceheader> <sourceframes> <outputtrajectory>"
	call getarg(1,headerfile)
	call getarg(2,framesfile)
	call getarg(3,outfile)

	write(0,"(A,A)") "Reading header from file : ",headerfile
	call openhis(headerfile,15)
        success = readheader()
	if (success.EQ.-1) stop "Couldn't read header from first file!"
	write(0,*) "Writing header of new file..."
	if (writeheader(outfile,14,0).EQ.-1) stop "Failed to write new history file header!"
	close(15)
	
	write(0,"(A,A)") "Reading frames from file : ",framesfile
	call openhis(framesfile,15)
	nframes = 0
10	success = readframe()
	if (success.NE.0) goto 15
	nframes = nframes + 1
	if (MOD(nframes,100).EQ.0) write(0,*) nframes
	! Write out the frame...
	success = writeframe()
	 if (success.EQ.-1) stop "Failed to write frame!"
	goto 10
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Read/write ",nframes,"from file."
	  
20      close(14)
	stop "Finished."

	end program addheader
