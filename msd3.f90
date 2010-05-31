!	** msd **
!	Calculate the mean square displacement of molecules (COM)

	program msd3
	use dlprw; use utility
	implicit none
	character*80 :: hisfiles(10),outfile,basename,resfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,m,s1,m1,nframes,success,nargs,baselen,file
	integer :: maxlength, count1, totmols, maxstats, framecount
	integer :: iargc, nhisfiles
	logical :: nomoreframes
	real*8, allocatable :: msd(:,:), msd_sp(:,:), x(:,:), y(:,:), z(:,:), dx(:,:), dy(:,:), dz(:,:)
	real*8, allocatable :: msd_stats(:,:,:)
	real*8, allocatable :: minimum(:,:),maximum(:,:),avg(:,:),sd(:,:),avgrij2(:,:)
	integer, allocatable :: msd_norm(:,:), msd_stats_norm(:)
	real*8 :: tx,ty,tz,rij2,deltat, total, num, deltax, deltay, deltaz

	nargs = iargc()
	if (nargs.lt.5) stop "Usage : msd3 <DLP OUTPUTfile> <delta t> <maxlength> <maxstats> <HISfile> [HISfile...]"
	call getarg(1,outfile)
	call getarg(2,temp); read(temp,"(F10.6)") deltat
	call getarg(3,temp); read(temp,"(I6)") maxlength
	call getarg(4,temp); read(temp,"(I6)") maxstats

	do n=5,nargs
	  call getarg(n,hisfiles(n-4))
	end do
	nhisfiles = nargs-4

	! Open and check the files...
	if (outinfo(outfile,1).EQ.-1) goto 798

	totmols = 0
	do n=1,nspecies
	  totmols = totmols + s_nmols(n)
	end do

	allocate (msd(maxlength,totmols))
	allocate (msd_sp(maxlength,nspecies))
	allocate (msd_norm(maxlength,nspecies))
	allocate (x(maxlength,totmols))
	allocate (y(maxlength,totmols))
	allocate (z(maxlength,totmols))
	allocate (dx(maxlength,totmols))
	allocate (dy(maxlength,totmols))
	allocate (dz(maxlength,totmols))
	if (maxstats.gt.0) then
	  allocate(msd_stats(maxlength,nspecies,maxstats), avgrij2(maxlength,nspecies), msd_stats_norm(maxlength))
	  allocate(minimum(maxlength,nspecies),maximum(maxlength,nspecies),sd(maxlength,nspecies), avg(maxlength,nspecies))
	  msd_stats = 0.0
	  msd_stats_norm = 0
	end if

	msd = 0.0
	msd_norm = 0

	! Loop over files
	do file=1,nhisfiles

	  x = 0.0
	  y = 0.0
	  z = 0.0
	  dx = 0.0
	  dy = 0.0
	  dz = 0.0
	  framecount = 0

	  write(0,"(a,i4,a,a)") "Moving to history file ",file,":",hisfiles(file)
	  close(10)
	  call openhis(hisfiles(file),10)
	  if (readheader().EQ.-1) then
	    write(0,*) "Failed to read header. Skipping..."
	    cycle
	  end if

	  nomoreframes = .false.

	  ! First, populate array with 'maxlength' frames
	  do nframes=1,maxlength-1
	    success = readframe()
	    framecount = framecount + 1
	    if (mod(framecount,100).eq.0) write(0,*) framecount,"(buffering)"
	    if (success.ne.0) exit
	    call calc_com
	    ! Store the initial positions
	    count1 = 0
	    do s1=1,nspecies
	      do m1=1,s_nmols(s1)
	        count1 = count1 + 1
	        x(nframes,count1) = comx(s1,m1)
	        y(nframes,count1) = comy(s1,m1)
	        z(nframes,count1) = comz(s1,m1)
		! Determine delta
		if (nframes.gt.1) then
		  call pbc(comx(s1,m1),comy(s1,m1),comz(s1,m1),x(nframes-1,count1),y(nframes-1,count1),z(nframes-1,count1),tx,ty,tz)
		  dx(nframes,count1) = x(nframes-1,count1) - tx
		  dy(nframes,count1) = y(nframes-1,count1) - ty
		  dz(nframes,count1) = z(nframes-1,count1) - tz
		end if
	      end do !m1
	    end do !s1
	  end do

	  ! Continue accumulating data and reading new frames from file until we get EOF or error

	  do
	    ! Read a new frame....
	    if (nomoreframes) then
	      if (nframes.eq.1) exit
	    else
	      success = readframe()
	      if (success.ne.0) then
	        ! No more frames to read, but continue to shuffle data until we have nothing left to work on
	        nomoreframes = .true.
		write(0,*) "Emptying buffer..."
	      end if
	    end if

	    if (.not.nomoreframes) then
	      framecount = framecount + 1
	      if (mod(framecount,100).eq.0) write(0,*) framecount

	      ! ... and store new frame data
	      call calc_com()
	      count1 = 0
	      do s1=1,nspecies
		do m1=1,s_nmols(s1)
		  count1 = count1 + 1
		  x(nframes,count1) = comx(s1,m1)
		  y(nframes,count1) = comy(s1,m1)
		  z(nframes,count1) = comz(s1,m1)
		  ! Determine delta
		  call pbc(comx(s1,m1),comy(s1,m1),comz(s1,m1),x(nframes-1,count1),y(nframes-1,count1),z(nframes-1,count1),tx,ty,tz)
		  dx(nframes,count1) = x(nframes-1,count1) - tx
		  dy(nframes,count1) = y(nframes-1,count1) - ty
		  dz(nframes,count1) = z(nframes-1,count1) - tz
		end do !m1
	      end do !s1
	    else
	      nframes = nframes - 1
	    end if


	    ! Accumulate data for stored positions
	    avgrij2 = 0.0
	    count1 = 0
	    do s1=1,nspecies
	      do m1=1,s_nmols(s1)
		count1 = count1 + 1
		deltax=0.0
		deltay=0.0
		deltaz=0.0

		do n=2,nframes
		  deltax = deltax + dx(n,count1)
		  deltay = deltay + dy(n,count1)
		  deltaz = deltaz + dz(n,count1)

		  ! Calculate r**2 and accumulate average (for stats)
		  rij2=(deltax*deltax + deltay*deltay + deltaz*deltaz)
		  avgrij2(n-1,s1) = avgrij2(n-1,s1) + rij2

		  ! Increment msd array
		  msd(n-1,count1) = msd(n-1,count1) + rij2
		  msd_norm(n-1,s1) = msd_norm(n-1,s1) + 1

		end do
	      end do
	    end do

	    ! ...accumulate statistics?....
	    if (1.le.maxstats) then
	      do n=2,nframes
		msd_stats_norm(n-1) = msd_stats_norm(n-1) + 1
		do s1=1,nspecies
		  msd_stats(n-1,s1,msd_stats_norm(n-1)) = avgrij2(n-1,s1)
		end do
	      end do
	    end if

	    ! ...and shuffle position data back one box in the arrays
	    do n=2,maxlength
	      x(n-1,:) = x(n,:)
	      y(n-1,:) = y(n,:)
	      z(n-1,:) = z(n,:)
	      dx(n-1,:) = dx(n,:)
	      dy(n-1,:) = dy(n,:)
	      dz(n-1,:) = dz(n,:)
	    end do

	  end do 

	end do

	goto 800

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
800	write(0,*) "Finished."
801	write(0,*) ""

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (outfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="rdfresults."
	  baselen=11
	else
	  basename=outfile(1:baselen)
	endif

	! Accumulate molecular MSDs into species MSDs, normalise, work out stats
	msd_sp = 0.0
	do n=1,maxlength-1
	  count1 = 0
	  do s1=1,nspecies

	    ! Sum molecular MSDs
	    do m1=1,s_nmols(s1)
	      count1 = count1 + 1
	      msd_sp(n,s1) = msd_sp(n,s1) + msd(n,count1)
	    end do

	    ! Normalise
	    msd_sp(n,s1) = msd_sp(n,s1) / real(msd_norm(n,s1))

	    ! Calculate statistics
	    if (maxstats.gt.0) then
	      ! Normalise
	      num = msd_stats_norm(n)
	      msd_stats(n,s1,1:num) = msd_stats(n,s1,1:num) / s_nmols(s1)
	      minimum(n,s1) = minval(msd_stats(n,s1,1:num))
	      maximum(n,s1) = maxval(msd_stats(n,s1,1:num))
	      avg(n,s1) = sum(msd_stats(n,s1,1:num)) / real(num)
	      ! Now for S.D.
	      total = 0.0d0
	      do m=1,num
		total = total + (msd_stats(n,s1,m) - avg(n,s1))**2
	      end do
	      sd(n,s1) = SQRT( total / real(num) )
	    end if

	  end do

	end do

	! Write data
	do s1=1,nspecies
	  open(unit=9,file=basename(1:baselen)//"msd"//char(48+s1),form="formatted",status="replace")
	  write(9,"(a)") "# DeltaT            MSD             SD             Min              Max             Avg"
	  do n=1,maxlength-1
	    if (maxstats.gt.0) then
	      write(9,"(8f16.9)") n*deltat,msd_sp(n,s1),sd(n,s1),minimum(n,s1),maximum(n,s1),avg(n,s1),msd_stats_norm(n)*1.0
	    else
	      write(9,"(8f16.9)") n*deltat,msd_sp(n,s1)
	    end if
	  end do
	  close(9)
	end do

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program msd3

