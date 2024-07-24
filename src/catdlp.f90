	! ###############
	! dlp concatenator
	! ###############

	program catdlp
	use dlprw
	implicit none
	character*80 :: files(100),outfile,arg,headerfile
	integer :: nargs,n,success,nframes,totframes
	integer :: iargc, nfiles = 0
	logical :: failed_header, skipfirst = .FALSE., skipfirstfirst = .FALSE., skip, altheader = .FALSE.

	nargs = iargc()
	if (nargs.lt.2) then
	  write(0,*) " Usage: catdlp <file1> <file2> .. <fileN> <target>"
	  write(0,*) "          -discardfirst  -discardotherfirst -header <file>"
	  stop
	end if
	n = 1
	do 
	  call getarg(n,arg)
	  select case (arg)
	    case ("-discardfirst")
	      skipfirst = .TRUE.
	      skipfirstfirst = .TRUE.
	    case ("-discardotherfirst")
	      skipfirst = .TRUE.
	      skipfirstfirst = .FALSE.
	    case ("-header")
	      altheader = .TRUE.
	      n = n + 1
	      call getarg(n,headerfile)
	    case default
	      nfiles = nfiles + 1
	      files(nfiles) = arg
	      write(6,"(a,i3,a,a80)") "File ", nfiles, " is : ", files(nfiles)
	  end select
	  if (n.eq.nargs) exit
	  n = n + 1
	end do
	! Last specified filename is output filename
	outfile = files(nfiles)
	nfiles = nfiles - 1
	if (nfiles.le.0) stop "Error: No input files were given." 
	write(6,*) "Output file will be : ", outfile
	! Test open the output file....
	open(unit=20, file=outfile, status='old', err=5)
	stop "Will not overwrite existing output file."

5	close(20)
	write(6,*) "Concatenating..."

	totframes = 0
	do n=1,nfiles
	  nframes = 0
	  failed_header = .FALSE.
	  ! Open the file and read in the file header (not used)
	  write(6,"(A,I2,A,A)") "Read/Write file ",n," : ",files(n)
	  call openhis(files(n),15)
          success = readheader()
	  if (n.EQ.1) then
	    ! First file, so write the header of the new file
	    if (success.EQ.-1) then
	      if (.not.altheader) stop "Couldn't read header of first file!"
	      call openhis(headerfile,15)
	      if (readheader().eq.-1) stop "Failed to read header from first file *and* alternate file"
	      call openhis(files(n),15)
	    end if
	    write(6,*) "Writing header of new file..."
	    if (writeheader(outfile,14,0).EQ.-1) stop "Failed to write new history file header!"
	  else
	    ! Its the nth file, so we don't really care if we can't read the header since it might
	    ! be a restarted trajectory. So we'll carry on anyway. However, set a flag here to check
	    ! when we read the first frame - if this also fails, then it's probably not a history file.
	    if (success.NE.0) then
	      write(6,*) "No header found in this file - restarted trajectory?"
	      rewind(dlpun_his)
	      failed_header = .TRUE.
	    else
	      write(6,"(A,I2)") "Read header of file ",n
	    end if
	  end if

	  ! Skip first frame of file?
	  skip = .FALSE.
	  if (skipfirst) then
	    ! Skip first frame of *first* file?
	    if ((n.ne.1) .or. ((n.eq.1).and.skipfirstfirst)) skip = .TRUE.
	  end if
	  if (skip) write(6,"(a,i3)") "Skipping first frame of file ", n

10	  success = readframe()
	  if ((nframes.EQ.0).AND.(failed_header).AND.(success.NE.0)) then
	    write(6,"(A,A)") "Failed to read header *and* first frame of file ",files(n)
	    stop "Will quit."
	  end if
	  if (success.NE.0) goto 15
	  if (skip) then
	    skip = .FALSE.
	    goto 10
	  end if

	  ! Increase counters
	  nframes = nframes + 1
	  totframes = totframes + 1

	  if (MOD(nframes,100).EQ.0) write(6,*) nframes
	  ! Write out the frame...
	  success = writeframe()
	  if (success.EQ.-1) stop "Failed to write frame!"
	  goto 10
15	  close(15)  ! end of file..
	  write(6,"(A,I5,A,A)") "Read/wrote ",nframes," from file ",files(n)
	end do
	  
20      close(14)
	write(6,"(a,i6,a)") "Read/wrote ",totframes," in total."
	stop "Finished concatenation."

	end program catdlp
