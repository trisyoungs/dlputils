	! ###############
	! dlp chunker
	! ###############

	program chunkdlp
	use dlprw
	implicit none
	character*80 :: hisfile,basename,outfile,framefile,temp
	integer :: nargs,n,i,m,o,success,nframes,nfiltered,tth,th,hun,ten,units
	integer :: iargc, framefirst, framelast, chunksize, nchunks, baselen

	nargs = iargc()
	if ((nargs.ne.3).and.(nargs.ne.5)) stop "Usage : chunkdlp <input HISfile> <basename (will be appended by '_XXX.HISu'> <startframe> <chunksize> <endframe>"
	call getarg(1,hisfile)
	call getarg(2,basename)
	call getarg(3,temp); read(temp,"(i7)") framefirst
	call getarg(4,temp); read(temp,"(i7)") chunksize
	call getarg(5,temp); read(temp,"(i7)") framelast
	write(6,"(a,i6,a,i6,a,i6)") "Using frame range ",framefirst," to ",framelast," with a chunksize of ", chunksize

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (basename(n:n).ne." ") then
	    baselen=n
	    goto 52
	  endif
	end do
52     if (baselen.EQ.-1) then
	  basename="hischunk"
	  baselen=8
	endif

	write(6,*) "Chunking..."
	nchunks = 0
	nframes = 0
	! Open the file and read in the file header (not used)
	call openhis(hisfile,15)
        success = readheader()
	if (success.EQ.-1) stop "Couldn't read header of first file!"

	write(6,*) "Moving to frame ", framefirst
	nfiltered = 0
10	success = readframe()
	if (success.NE.0) then
	  write(0,"(A,i6,A)") "Failed to read frame ", nframes," from file"
	  stop "Quit."
	end if
	nframes = nframes + 1
	if (nframes.lt.framefirst) goto 10
	if (nfiltered.eq.0) then
	  write(6,*) "Beginning chunks from frame ",nframes
	  ! Construct outfile name
	  tth = nframes / 10000; n = nframes - tth*10000
	  th = n / 1000; n = n - th*1000
	  hun = n / 100; n = n - hun*100
	  ten = n / 10; n = n - ten*10
	  units = n
	  outfile = basename(1:baselen)//"_"//CHAR(48+tth)//CHAR(48+th)//CHAR(48+hun)//CHAR(48+ten)//CHAR(48+units)//".HISu"
	write (0,*) outfile
	  write(6,*) "Writing header of new file "
	  if (writeheader(outfile,14,0).EQ.-1) stop "Failed to write new history file header!"
	  nchunks = nchunks + 1
	end if

	! Write out the frame...
	success = writeframe()
	if (success.EQ.-1) stop "Failed to write frame!"
	nfiltered = nfiltered + 1

	if (nfiltered.eq.chunksize) then
	  nfiltered = 0
	  close(14)
	end if

	if (nframes.eq.framelast) goto 15
	goto 10

14	write(0,*) "Found end of frames file."
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Filtered ",nfiltered," frames from file."
	goto 20

17	stop "Error reading from frames file."
	  
20      close(14)
	stop "Finished chunking."

	end program chunkdlp
