!	** msd **
!	Calculate the mean square displacement of molecules (COM)

	program msdprog
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,m,s1,m1,nframes,success,nargs,baselen
	integer :: framestodo, count1, totmols, framestodiscard = 0
	integer :: iargc
	real*8, allocatable :: lastx(:), lasty(:), lastz(:), msd(:,:)
	real*8, allocatable :: msd_x(:,:), msd_y(:,:), msd_z(:,:), msd_sp(:,:)
	real*8, allocatable :: dx(:), dy(:), dz(:)
	integer, allocatable :: orn(:,:)
	real*8 :: tx,ty,tz,rij2,deltat

	nargs = iargc()
	if (nargs.lt.4) stop "Usage : msd <HISTORYfile> <OUTPUTfile> <delta t> <framestodo> [framestodiscard]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(F10.6)") deltat
	call getarg(4,temp); read(temp,"(I6)") framestodo
	if (nargs.gt.4) then
	  call getarg(5,temp); read(temp,"(I6)") framestodiscard
	end if

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
	! Now, read in the history header so that we have cell()
	if (readheader().EQ.-1) goto 799

	totmols = 0
	do n=1,nspecies
	  totmols = totmols + s_nmols(n)
	end do

	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....

	call calc_com

	allocate (orn(framestodo,nspecies))
	allocate (msd(framestodo,totmols))
	allocate (msd_sp(framestodo,nspecies))
	allocate (msd_x(framestodo,totmols))
	allocate (msd_y(framestodo,totmols))
	allocate (msd_z(framestodo,totmols))
	allocate (lastx(totmols))
	allocate (lasty(totmols))
	allocate (lastz(totmols))
	allocate (dx(totmols))
	allocate (dy(totmols))
	allocate (dz(totmols))

	orn = 0
	msd = 0.0
	msd_sp = 0.0
	msd_x = 0.0
	msd_y = 0.0
	msd_z = 0.0
	
	! Store the initial positions
	count1 = 0
	do s1=1,nspecies
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1

	    lastx(count1) = comx(s1,m1)
	    lasty(count1) = comy(s1,m1)
	    lastz(count1) = comz(s1,m1)

	  end do !m1
	end do !s1

100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	if (framestodiscard.ne.0) then
	  framestodiscard = framestodiscard - 1
	  goto 101
	end if
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	call calc_com

	! Calculate displacements of species between this frame and the last
	count1 = 0
	do s1=1,nspecies
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1

	    call pbc(comx(s1,m1),comy(s1,m1),comz(s1,m1),lastx(count1),lasty(count1),lastz(count1),tx,ty,tz)
	    dx(count1) = lastx(count1) - tx
	    dy(count1) = lasty(count1) - ty
	    dz(count1) = lastz(count1) - tz

	  end do !m1
	end do !s1

	! Store current positions for use in next iteration
	count1 = 0
	do s1=1,nspecies
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1

	    lastx(count1) = comx(s1,m1)
	    lasty(count1) = comy(s1,m1)
	    lastz(count1) = comz(s1,m1)

	  end do !m1
	end do !s1


	! Update averages over origins
	do n=1,nframes
	  count1 = 0
	  do s1=1,nspecies
	    do m1=1,s_nmols(s1)
	      count1 = count1 + 1

	      ! Origin COM positions given by rxyz(m,count1)
	      ! 'To' COM position given by ryxz(nframes,count1)

	      ! Increase total msd_xyz
	      msd_x(n,count1) = msd_x(n,count1) + dx(count1)
	      msd_y(n,count1) = msd_y(n,count1) + dy(count1)
	      msd_z(n,count1) = msd_z(n,count1) + dz(count1)

	      ! Calculate r**2
	      rij2=(msd_x(n,count1))**2 + (msd_y(n,count1))**2 + (msd_z(n,count1))**2

	      ! Increment msd array
	      msd(nframes-(n-1),count1) = msd(nframes-(n-1),count1) + rij2
	      orn(nframes-(n-1),s1) = orn(nframes-(n-1),s1) + 1

	    end do !m1
	  end do !s1
	end do !frames


	if (nframes.EQ.framestodo) goto 800
	! Next frame
	goto 101

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended before framestodo was fulfilled..."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",framestodo,")"
	goto 801
800	write(0,*) "Framestodo was fulfilled."
801	write(0,*) ""

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="rdfresults."
	  baselen=11
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=9,file=basename(1:baselen)//"msd",form="formatted",status="replace")
	do n=1,nframes
	  count1 = 0
	  do s1=1,nspecies
	    do m1=1,s_nmols(s1)
	      count1 = count1 + 1

	      ! Sum over number of molecules
	      msd_sp(n,s1) = msd_sp(n,s1) + msd(n,count1)

	    end do

	    ! Average over number of molecules and origins
	    msd_sp(n,s1) = msd_sp(n,s1) / real(orn(n,s1))
	    !msd_sp(n,s1) = msd_sp(n,s1) / s_nmols(s1)
	
	  end do

	  ! Write the data...
	  write(9,"(8f15.9)") n*deltat,(msd_sp(n,m),m=1,nspecies)
	end do
	close(9)

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program msdprog

