	! ###############
	! dlp concatenator
	! ###############

	program catdlp
	use dlprw
	implicit none
	character*80 :: files(100),outfile,arg
	integer :: nargs,n,success,nframes,totframes
	integer :: iargc, nfiles = 0
	logical :: failed_header, skipfirst = .FALSE., skipfirstfirst = .FALSE.

	nargs = iargc()
	if (nargs.lt.2) then
	  write(0,*) " Usage: catdlp <file1> <file2> .. <fileN> <target>"
	  write(0,*) "          -skipfirst  -skipotherfirst"
	end if
	do n=1,nargs
	  call getarg(n,arg)
	  select case (arg)
	    case ("-skipfirst")
	      skipfirst = .TRUE.
	      skipfirstfirst = .TRUE.
	    case ("-skipotherfirst")
	      skipfirst = .TRUE.
	      skipfirstfirst = .FALSE.
	    case default
	      nfiles = nfiles + 1
	      files(nfiles) = arg
	      write(6,"(a,i3,a,a80)") "File ", nfiles, " is : ", files(nfiles)
	  end select
	end do
	! Last specified filename is output filename
	outfile = files(nfiles)
	nfiles = nfiles - 1
	if (nfiles.le.0) stop "Error: No input files were given." 
	write(6,*) "Output file will be : ", outfile

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
	    if (success.EQ.-1) stop "Couldn't read header of first file!"
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

10	  success = readframe()
	  if ((nframes.EQ.0).AND.(failed_header).AND.(success.NE.0)) then
	    write(6,"(A,A)") "Failed to read header *and* first frame of file ",files(n)
	    stop "Will quit."
	  end if
	  if (success.NE.0) goto 15

	  ! Skip first frame?
	  if ((nframes.eq.0).and.skipfirst) then
	    ! Skip first frame of *first* file?
	    if ((n.ne.1) .or. ((n.eq.1).and.skipfirstfirst)) then
	      write(6,"(a,i3)") "Skipping first frame of file ", n
	      continue
	    end if
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