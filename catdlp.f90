	! ###############
	! dlp concatenator
	! ###############

	program catdlp
	use dlprw
	implicit none
	character*80 :: files(100),outfile
	integer :: nargs,n,i,m,o,success,nframes,totframes
	integer :: iargc
	logical :: failed_header

	nargs = iargc()
	if (nargs.LT.2) stop " Not enough arguments supplied: <file1> <file2> .. <fileN> <target>"
	call getarg(nargs,outfile)
	do n=1,nargs-1
	  call getarg(n,files(n))
	end do

	write(6,*) "Concatenating..."

	totframes = 0
	do n=1,nargs-1
	  nframes = 0
	  failed_header = .FALSE.
	  ! Open the file and read in the file header (not used)
	  write(6,"(A,I2,A,A)") "Read/Write file ",n," : ",files(n)
	  call openhis(files(n),15)
          success = readheader()
	  if (n.EQ.1) then
	    ! First file, so write the header of the new file
	    if (success.EQ.-1) stop "Couldn't read header of first file!"
	    write(0,*) "Writing header of new file..."
	    if (writeheader(outfile,14,0).EQ.-1) stop "Failed to write new history file header!"
	  else
	    ! Its the nth file, so we don't really care if we can't read the header since it might
	    ! be a restarted trajectory. So we'll carry on anyway. However, set a flag here to check
	    ! when we read the first frame - if this also fails, then it's probably not a history file.
	    if (success.NE.0) then
	      write(0,*) "No header found in this file - restarted trajectory?"
	      rewind(dlpun_his)
	      failed_header = .TRUE.
	    else
	      write(0,"(A,I2)") "Read header of file ",n
	    end if
	  end if

10	  success = readframe()
	  if ((nframes.EQ.0).AND.(failed_header).AND.(success.NE.0)) then
	    write(0,"(A,A)") "Failed to read header *and* first frame of file ",files(n)
	    stop "Will quit."
	  end if
	  if (success.NE.0) goto 15
	  nframes = nframes + 1
	  totframes = totframes + 1
	  if (MOD(nframes,100).EQ.0) write(0,*) nframes
	  ! Write out the frame...
	  success = writeframe()
	  if (success.EQ.-1) stop "Failed to write frame!"
	  goto 10
15	  close(15)  ! end of file..
	  write(0,"(A,I5,A,A)") "Read/wrote ",nframes," from file ",files(n)
	end do
	  
20      close(14)
	write(0,"(a,i6,a)") "Read/wrote ",totframes," in total."
	stop "Finished concatenation."

	end program catdlp
