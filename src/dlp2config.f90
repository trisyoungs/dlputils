	! ###############
	! dlp2config
	! ###############

	program dlp2config
	use dlprw
	implicit none
	character*80 :: hisfile,temp,altheaderfile,destfile
	integer :: nargs,n,i,m,success,nframes,targetframe
	integer :: iargc,tth,th,hun,ten,units
	logical :: altheader = .false.

	nargs = iargc()
	if (nargs.lt.2) stop "Usage: dlp2config <HISTORYfile> <frame no. | -1 for all or every> [altheader]"
	call getarg(1,hisfile)
	call getarg(2,temp)
	read(temp,"(I7)") targetframe
	if (nargs.eq.3) then
	  call getarg(3,altheaderfile)
	  altheader = .true.
	end if

	! Open the file and read in the file header
	call openhis(hisfile,15)
        success = readheader()
	if (success.EQ.-1) then
	  write(0,*) "Couldn't read header of first file!"
	  write(0,*) "Restarted trajectory?"
	  rewind(dlpun_his)
          if (altheader) then
            close(dlpun_his)
            call openhis(altheaderfile,15)
            if (readheader().EQ.-1) stop "Couldn't read header from alternative file."
            close(dlpun_his)
            call openhis(hisfile,15)
	    write(0,*) "Successfully got header information from alternative file."
	  else
	    stop "Provide the name of an alternate HIS file (with a suitable header) on the command line."
	  end if
	end if

	write(6,*) "Seeking..."

	nframes = 0
10	success = readframe()
	if (success.NE.0) goto 15
	nframes = nframes + 1
	if (MOD(nframes,100).EQ.0) write(0,*) nframes
	if ((nframes.ne.targetframe).and.(targetframe.gt.0)) goto 10
	if ((targetframe.lt.0).and.(mod(nframes,abs(targetframe)).ne.0)) goto 10

	! Copy history file configuration to CONFIG arrays
	do n=1,natms
	  cfgxyz(n,1) = xpos(n)
	  cfgxyz(n,2) = ypos(n)
	  cfgxyz(n,3) = zpos(n)
	  cfgvel(n,1) = xvel(n)
	  cfgvel(n,2) = yvel(n)
	  cfgvel(n,3) = zvel(n)
	  cfgfrc(n,1) = xfor(n)
	  cfgfrc(n,2) = yfor(n)
	  cfgfrc(n,3) = zfor(n)
	  cfgname(n) = atmname(n)
	end do
	cfgcell = cell
	cfgheader = dlpheader
	
	! Write out the configuration
	write(0,*) nframes
	! Construct destfile name
	tth = nframes / 10000; n = nframes - tth*10000
	th = n / 1000; n = n - th*1000
	hun = n / 100; n = n - hun*100
	ten = n / 10; n = n - ten*10
	units = n
	
	destfile = "frame"//CHAR(48+tth)//CHAR(48+th)//CHAR(48+hun)//CHAR(48+ten)//CHAR(48+units)//".CONFIG"
	write(0,*) "Frame data will be written to :",destfile

	call writeconfig(destfile)

	if (targetframe.lt.0) goto 10

15	close(15)  ! end of file..
	  
	stop "Finished."

	end program dlp2config
