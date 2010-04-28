	! ###########
	! dlp thinner
	! ###########

	program thindlp
	use dlprw
	implicit none
	character*80 :: infile, outfile, temp
	integer :: nargs,n,i,m,o,success,nframes,nwritten
	integer :: iargc, skipatstart, nskip
	logical :: failed_header

	nargs = iargc()
	if (nargs.ne.4) stop "Usage is : dlpthin <source.HISu> <nSkipAtStart> <frameSkip> <dest.HISu>"
	call getarg(nargs,outfile)
	call getarg(1,infile)
	call getarg(4,outfile)
	call getarg(2,temp)
	read(temp,"(i5)") skipatstart
	call getarg(3,temp)
	read(temp,"(i5)") nskip
	

	write(6,*) "Thinning..."

	nframes = 0
	nwritten = 0
	failed_header = .FALSE.
	! Open the file and read in the file header (not used)
	write(6,"(A,A)") "Opening source file ",infile
	call openhis(infile,15)
        success = readheader()
	! First file, so write the header of the new file
	if (success.EQ.-1) then
	  write(0,*) "Couldn't read header of input file - no header will be written to output file..."
	  rewind(dlpun_his)
	else
	  write(0,*) "Writing header of new file..."
	  if (writeheader(outfile,14,0).EQ.-1) stop "Failed to write new history file header!"
	end if

10	success = readframe()
	if ((nframes.EQ.0).AND.(failed_header).AND.(success.NE.0)) then
	  write(0,"(A,A)") "Failed to read header *and* first frame of file ",infile
	  stop "Will quit."
	end if
	if (success.NE.0) goto 15
	nframes = nframes + 1
	if (MOD(nframes,100).EQ.0) write(0,*) nframes
	if (nframes.lt.skipatstart) goto 10
	! Write out the frame...
	if (mod((nframes - skipatstart),nskip).eq.0) then
	  success = writeframe()
	  if (success.EQ.-1) stop "Failed to write frame!"
	  nwritten = nwritten + 1
	end if
	goto 10
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Read ",nframes," from ",infile
	write(0,"(A,I5,A,A)") "Wrote ",nwritten," to ",outfile
	  
20      close(14)
	stop "Finished thinning."

	end program thindlp
