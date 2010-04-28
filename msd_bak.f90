!	** msd **
!	Calculate the mean square displacement of molecules (COM)

	program msdprog
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,m,s1,m1,nframes,success,nargs
	integer :: framestodo, count1, totmols
	integer :: iargc
	real*8, allocatable :: rx(:,:), ry(:,:), rz(:,:), msd(:,:)
	integer, allocatable :: orn(:,:)
	real*8 :: tx,ty,tz,rij2,deltat

	nargs = iargc()
	if (nargs.NE.4) then
	  write(0,"(A)") "Usage : msd <DLP HISTORYfile> <DLP OUTPUTfile> <delta t> <framestodo>"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(F10.6)") deltat
	call getarg(4,temp); read(temp,"(I6)") framestodo

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
	if (success.EQ.-1) goto 799  ! File error....

	call calc_com

	allocate (orn(framestodo,nspecies))
	allocate (msd(framestodo,nspecies))
	allocate (rx(0:framestodo,totmols))
	allocate (ry(0:framestodo,totmols))
	allocate (rz(0:framestodo,totmols))

	orn = 0
	msd = 0.0
	
	! Store the initial COM positions
	count1 = 0
	do s1=1,nspecies
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1

	    rx(0,count1) = comx(s1,m1)
	    ry(0,count1) = comy(s1,m1)
	    rz(0,count1) = comz(s1,m1)

	  end do !m1
	end do !s1

100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.EQ.-1) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	call calc_com

	! Store current configuration COM positions
	count1 = 0
	do s1=1,nspecies
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1

	    rx(nframes,count1) = comx(s1,m1)
	    ry(nframes,count1) = comy(s1,m1)
	    rz(nframes,count1) = comz(s1,m1)

	  end do !m1
	end do !s1

	! Update averages over origins
	do n=0,nframes-1
	  count1 = 0
	  do s1=1,nspecies
	    do m1=1,s_nmols(s1)
	      count1 = count1 + 1

	      ! Origin COM positions given by rxyz(m,count1)
	      ! 'To' COM position given by ryxz(nframes,count1)

	      call pbc(rx(n,count1),ry(n,count1),rz(n,count1), &
	      &    rx(nframes,count1),ry(nframes,count1),rz(nframes,count1),tx,ty,tz)
	      rij2=(tx-rx(nframes,count1))**2 + (ty-ry(nframes,count1))**2 + (tz-rz(nframes,count1))**2
	      ! Add the distance to the msd array
	      !write(0,"(A,2F8.4)") "rij2 =",rij2
	      msd(nframes - n,s1) = msd(nframes - n,s1) + rij2
	      orn(nframes - n,s1) = orn(nframes - n,s1) + 1

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

	! Perform averaging over origins and molecules
	open(unit=9,file="msd.dat",form="formatted",status="replace")
	do n=1,nframes
	  do s1=1,nspecies

	    ! Average over number of molecules
	    msd(n,s1) = msd(n,s1) / (s_nmols(s1) * (nframes - (n-1)))

	  end do
	  ! Write the data...
	  write(9,"(5F10.4)") n*deltat,(msd(n,m),m=1,nspecies)
	end do
	close(9)

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program msdprog

