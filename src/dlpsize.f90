	! ###############
	! dlp history sizer
	! ###############

	program dlpsize
	use dlprw
	implicit none
	character*80 :: hisfile
	integer :: nargs,n,i,m,o,success,nframes
	integer :: iargc
	logical :: failed_header

	nargs = iargc()
	if (nargs.LT.1) stop "Usage: dlpsize <HISTORYfile>"
	call getarg(1,hisfile)

	write(6,*) "Counting..."

	nframes = 0
	failed_header = .FALSE.
	! Open the file and read in the file header (not used)
	write(6,"(A,I2,A,A)") "Read/Write file ",n," : ",hisfile
	call openhis(hisfile,15)
        success = readheader()
	if (success.NE.0) then
	  write(0,*) "Success = ",success
	  write(0,*) "No header found in this file - restarted trajectory?"
	  rewind(dlpun_his)
	end if

10	success = readframe()
	if ((nframes.EQ.0).AND.(success.NE.0)) then
	  stop "Failed to read first frame of file."
	end if
	if (success.NE.0) goto 15
	nframes = nframes + 1
	!if (MOD(nframes,100).EQ.0) write(0,*) nframes
	goto 10
15	close(15)  ! end of file..
	write(0,"(A,I10,A,A)") "Successfully read ",nframes," from file."
	  
	stop "Finished count."

	end program dlpsize
