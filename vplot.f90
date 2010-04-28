	program vplot
	implicit none
	integer :: PGOPEN, nargs, n, m, nbins, nframes, x11id, i, currentframe
	character*80 :: resfile,fffile,command
	real, allocatable :: hist(:,:), binx(:)
	real*8 :: binwidth, ymax, xmax, temp, storedymax

	nargs = iargc()
	if (nargs.NE.1) stop "Usage: vplot <results file>"
	! Get command line args
	call getarg(1,resfile)

	! Read in histogram data from file.
	! Determine number of 'frames'
	open(unit=11,file=resfile,form="formatted",status="old")
	nframes = 0
10	read (11,*,end=15,err=15) nbins,binwidth
	nframes = nframes + 1
	do n=1,nbins
	  read(11,*) temp
	end do
	goto 10
15	write(0,*) "Found ",nframes," histogram sets in file."
	write(0,*) "NBins = ",nbins
	write(0,*) "Binwidth = ",binwidth
	allocate(hist(nframes,nbins))
	allocate(binx(nbins))
	storedymax = 0.0
	rewind(11)
	do n=1,nframes
	  read(11,*) nbins,binwidth
	  do m=1,nbins
	    read(11,"(i)") i
	    hist(n,m) = float(i)
	    if (storedymax.lt.hist(n,m)) storedymax = hist(n,m)
	  end do
	end do
	write(0,*) "Read in all data successfully."
	! Calculate plot limits the binx array
	do n=1,nbins
	  binx(n) = binwidth*(real(n)-0.5)
	end do
	xmax = binwidth * nbins
	ymax = storedymax
	  
	! Initialise PGPLOT
        ! Open the default (X11) device
        x11id = PGOPEN('/XWIN')
        if (x11id.LE.0) stop "Failed to initialise PGPLOT."
        call PGASK(.false.)
        call PGSCR(0,1.0,1.0,1.0)       ! Set Background colour to white
        call PGSCR(1,0.0,0.0,0.0)       ! Set main pen to black

	! Draw first frame on screen
	call draw(nbins,xmax,ymax,binx,hist(1,:))
	currentframe = 1

	! The main program loop is here - read commands and process them
100	call draw(nbins,xmax,ymax,binx,hist(currentframe,:))
	write(0,*) "Current frame ",currentframe
	read(5,"(A)") command
	select case (command)
	  case ("run")
	    do n=1,nframes
	      write(0,*) "Frame ",n
	      call draw(nbins,xmax,ymax,binx,hist(n,:))
	      call mswait(500)
	    end do
	  case ("next","n")
	    currentframe = currentframe + 1
	    if (currentframe.gt.nframes) currentframe = 1
	  case ("prev","p")
	    if (currentframe.eq.0) currentframe = nframes
	    currentframe = currentframe - 1
	  case ("goto","g")
	    write(0,*) "Enter frame to go to : min = 1, max = ",nframes
	    read(5,"(A)") command
	    write(currentframe,"(i10)") command
	  case ("ymax")
	    write(0,*) "Enter new Y maximum:"
	    read(5,"(A)") command
	    write(ymax,"(f20.14)") command
	  case ("xmax")
	    write(0,*) "Enter new X maximum:"
	    read(5,"(A)") command
	    write(xmax,"(f20.14)") command
	  case ("reset")
	    ymax = storedymax
	    xmax = nbins * binwidth
	  case ("help","h")
	    write(0,*) "Commands available are:",""
	    write(0,*) "    run          Display all frames one after the other"
	    write(0,*) "    delay        Change the delay between frames",""
	    write(0,*) "    next or n    Skip to the next frame"
	    write(0,*) "    prev or p    Skip to the previous frame"
	    write(0,*) "    goto or g    Jump to a specific frame number"
	    write(0,*) "    ymax         Set a new Y maximum"
	    write(0,*) "    xmax         Set a new X maximum"
	    write(0,*) "    reset        Reset maxima to initial values"
	    write(0,*) "    quit or q    Exit"
	  case ("quit","q")
	    goto 999
	  case default
	    write(6,*) "Unrecognised command : ",command
	end select
	goto 100

999	call PGEND
	end program vplot

	subroutine mswait(ms)
	implicit none
	integer :: ms, clock_start, clock_now, clock_rate
	! Get tick rate of clock
	call system_clock(count_rate=clock_rate)
	! Start timer
	call system_clock(count=clock_start)
	! Waste loop - get current tick time
10	call system_clock(count=clock_now)
	if (ms.gt.(clock_now-clock_start)) goto 10
	!call system_clock()
	end subroutine mswait

        subroutine printpage(x11id,psname)
        implicit none
        character*(*) :: psname
        integer :: PGOPEN, psid, x11id
        ! Open the postscript device
        if (psname.EQ."") then
          psid = PGOPEN("/CPS")
        else
          psid = PGOPEN(psname//"/CPS")
        end if
        if (psid.LE.0) stop "Couldn't open print device!"
        call draw
	call PGCLOS(psid)
	call PGSLCT(x11id)
        end subroutine printpage

	subroutine draw(nbins,xmax,ymax,xdata,ydata)
	implicit none
        integer :: n,id, i, nxy, nbins
	real*8 :: xmax, ymax
	real :: xdata(nbins), ydata(nbins)
	call PGERAS
	call PGSLW(2)
	call PGSCH(1.0)				! Set default character height
	call PGSCI(1)   		      ! Set pen to black
	! Define scale of page
	call PGENV(0.0,xmax,0.0,ymax,0,2)
        ! Define panel area
	call PGBIN(nbins,xdata,ydata,.true.)
	end subroutine draw
