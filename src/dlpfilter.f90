	! ###############
	! dlp filter
	! ###############

	program dlpfilter
	use dlprw
	implicit none
	character*80 :: hisfile,outfile,framefile,temp
	integer :: nargs,n,i,m,o,success,nframes,seek,nfiltered
	integer :: iargc, framefirst, framelast, frameskip
	logical :: fromfile = .true.

	nargs = iargc()
	if ((nargs.ne.3).and.(nargs.ne.5)) then
	  write(0,*) "Usage : dlpfilter <input HISfile> <output HISfile> <frame file | first skip last>"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,outfile)
	if (nargs.eq.3) then
	  call getarg(3,framefile)
	  write(6,*) "Using frames listed in supplied file...."
	else
	  call getarg(3,temp); read(temp,"(i7)") framefirst
	  call getarg(4,temp); read(temp,"(i7)") frameskip
	  call getarg(5,temp); read(temp,"(i7)") framelast
	  write(6,"(a,i6,a,i6,a,i6)") "Using frame range ",framefirst," to ",framelast," with a skip of ", frameskip
	  fromfile = .false.
	end if

	write(6,*) "Filtering..."
	nfiltered = 0

	! Open the frames file and read in the first target (if using a file)
	if (fromfile) then
	  open(unit=19,file=framefile,form="formatted",status="old")
	  read(19,"(I10)") seek
	else
	  seek = framefirst
	end if

	nframes = 0
	! Open the file and read in the file header (not used)
	call openhis(hisfile,15)
        success = readheader()
	if (success.EQ.-1) stop "Couldn't read header of first file!"
	write(0,*) "Writing header of new file..."
	if (writeheader(outfile,14,0).EQ.-1) stop "Failed to write new history file header!"

	write(0,"(A,I6,A)") "Seeking frame ",seek,"..."

10	success = readframe()
	if (success.NE.0) then
	  write(0,"(A,i6,A)") "Failed to read frame ", nframes," from file"
	  stop "Quit."
	end if
	nframes = nframes + 1

	if (nframes.EQ.seek) then
	  ! Write out the frame...
	  success = writeframe()
	  if (success.EQ.-1) stop "Failed to write frame!"
	  nfiltered = nfiltered + 1

	  ! Read in the next target frame (if using a file)
	  if (fromfile) then
	    read(19,"(I10)",err=17,end=14) seek
	  else
	    seek = seek + frameskip + 1
	    if (seek.gt.framelast) goto 15
	  end if
	  write(0,"(A,I6,A)") "Seeking frame ",seek,"..."
	end if
	goto 10

14	write(0,*) "Found end of frames file."
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Filtered ",nfiltered," frames from file."
	goto 20

17	stop "Error reading from frames file."
	  
20      close(14)
	stop "Finished filtering."

	end program dlpfilter
