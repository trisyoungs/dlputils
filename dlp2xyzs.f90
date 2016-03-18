	! ###############
	! dlp2xyzs
	! ###############

	program dlp2xyzs
	use dlprw
	implicit none
	character*80 :: hisfile,outfile
	character*4 :: frameid
	integer :: nargs,n,success,nframes,dot,space, interval
	integer :: iargc
	logical :: failed_header

	nargs = iargc()
	if (nargs.NE.2) stop "Usage: dlp2xyzs <HISTORYfile> <interval>"
	call getarg(1,hisfile)
	call getarg(2,frameid)
	read (frameid,"(i4)") interval

	! Open the file and read in the file header
	call openhis(hisfile,15)
        success = readheader()
	if (success.EQ.-1) then
	  write(0,*) "Couldn't read header of first file!"
	  write(0,*) "Restarted trajectory?"
	end if

	write(0,*) "Skip interval is ",interval

        ! Make up output basename
        dot = 0; space = 0
        do n=1,80
          if (hisfile(n:n).eq.".") dot = n       ! Pos of last dot
          if ((hisfile(n:n).eq." ").and.(space.eq.0)) space = n  ! Pos of first space
         end do
        if (dot.EQ.0) dot=space

	nframes = 0
10	success = readframe()
	if (success.NE.0) goto 15
	nframes = nframes + 1
	if (MOD(nframes,100).EQ.0) write(0,*) nframes
	if (mod(nframes-1,interval).ne.0) goto 10

	! Write out xyz file of these coordinates
	write(frameid,"(i4)") nframes
	do n=1,3
	  if (frameid(n:n).eq." ") frameid(n:n) = "0"
	end do
        outfile = hisfile(1:dot-1)//"-frame"//frameid//".xyz"
	open(unit=10,file=outfile,form="formatted",status="replace")

	write(10,*) natms
	write(10,"(a,9(f9.4,1x))") "Cell ", cell
	do n=1,natms
	  write(10,"(a8,2x,3(f12.8,1x))")  atmname(n),xpos(n),ypos(n),zpos(n)
	end do
	close(10)

	goto 10

15	close(15)  ! end of file..
	  
	stop "Finished."

	end program dlp2xyzs
