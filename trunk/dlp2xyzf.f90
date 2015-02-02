	! ################
	! dlp2xyzf converter
	! Assumes that the trajectory file contains successive frames from the simulation
	! ################

	program dlp2xyzf
	use dlprw
	implicit none
	character*80 :: infile,temp,fmfile
	integer :: fileform,count,limit,n,m,i,nargs,success
	integer :: iargc
	real*8, allocatable :: sx(:) ,sy(:) ,sz(:)

	nargs = iargc()
	if (nargs.NE.2) stop "Syntax : dlp2xyzf <HISTORYfile> <nframes or 'all'>"
	call getarg(1,infile)
	write(0,"(A,A)") "Infile  : ",infile
	call getarg(2,temp)
	if (temp(1:3).EQ."all") then
	  limit = -10
	  write(0,"(A,A)") "Frames  : ","All"
	else
	  read(temp,"(I6)") limit
	  write(0,"(A,I6)") "Frames  : ",limit
	end if
	  
	call openhis(infile,15)
	if (readheader().EQ.-1) stop "Failed to read HISTORY file header!"

	allocate(sx(natms)); allocate(sy(natms)); allocate(sz(natms))

	write(0,*) "Converting....."
	count = 0
10	success = readframe()
	if (success.NE.0) then
	  if (success.EQ.1) write(0,*) "End of HISTORY file encountered."
	  if (success.EQ.-1) write(0,*) "HISTORY file ended prematurely!"
	  write(0,*) "FM frames written :",count
	  goto 20
	end if
	count = count + 1
	if (count.GT.1) then
	  write(6,*) natms
	  write(6,*) "Frame", count
	  ! Write out the new 'frame' before storing the new forces....
	  do n=1,natms
	    write(6,"(a8,2x,6e16.7)") atmname(n), sx(n),sy(n),sz(n),xfor(n),yfor(n),zfor(n)
	  end do
	end if
	! Store the current forces...
	sx = xpos
	sy = ypos
	sz = zpos
	if (MOD(count,50).EQ.0) write(0,*) "Count : ",count
	if (count-1.NE.limit) goto 10
20      close(14)
	stop "Finished conversion."

	end program dlp2xyzf
