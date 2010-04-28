!	** cagecor **
!	Calculate the cage correlation function for species in the system
!	(Rabani, Gezelter, Berne - J. Chem. Phys, 1997, 107, 6867)

	program cagecor
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,sp1,sp2,m1,m2,baselen,nframes,success,nargs,c, l_0t, n_i_out
	integer :: framestodo, count, sp1end, sp2end, totmols, s1, s2, count1, count2
	integer, allocatable :: l_t(:,:,:), l_tsq(:,:),ccorn(:)
	integer :: iargc
	real*8, allocatable :: ccor(:)
	real*8 :: c1x,c1y,c1z,tx,ty,tz,rij,cagecut

	nargs = iargc()
	if (nargs.LT.6) then
	  write(0,"(A)") "Usage : cagecor <DLP HISTORYfile> <DLP OUTPUTfile> <centresp> <secondsp> <cagecutoff> <c> [framestodo]"
	  write(0,"(A)") "	 Set <centresp> or <secondsp> to zero to specify 'any molecule'"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(I4)") sp1
	call getarg(4,temp); read(temp,"(I4)") sp2
	call getarg(5,temp); read(temp,"(F10.8)") cagecut
	call getarg(6,temp); read(temp,"(I6)") c
	if (nargs.EQ.7) then
	  call getarg(7,temp); read(temp,"(I6)") framestodo
	else
	  framestodo = -1
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
	if (sp1.EQ.0) then
	  sp1 = 1
	  sp1end = nspecies
	  count1 = totmols
	else
	  sp1end = sp1
	  count1 = s_nmols(sp1)
	endif
	if (sp2.EQ.0) then
	  sp2 = 1
	  sp2end = nspecies
	  count2 = totmols
	else
	  sp2end = sp2
	  count2 = s_nmols(sp2)
	endif

	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.EQ.-1) goto 799  ! File error....

	call calc_com

	! Construct the initial neighbour lists...
	allocate (l_t(0:framestodo,count1,count2))
	allocate (l_tsq(0:framestodo,count1))
	allocate (ccor(framestodo))
	allocate (ccorn(framestodo))

	ccor = 0.0
	ccorn = 0
	
	count1 = 0
	do s1=sp1,sp1end
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1
	    l_tsq(0,count1) = 0
	    c1x=comx(s1,m1)
	    c1y=comy(s1,m1)
	    c1z=comz(s1,m1)
	    count2 = 0
	    do s2=sp2,sp2end
	      do m2=1,s_nmols(s2)
		count2 = count2 + 1
		if ((s1.EQ.s2).AND.(m1.EQ.m2)) then
		  l_t(0,count1,count2) = 0
		else
		  call pbc(comx(s2,m2),comy(s2,m2),comz(s2,m2),c1x,c1y,c1z,tx,ty,tz)
		  rij=sqrt( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )
		  ! Heaviside function
		  if (rij.LT.cagecut) then
		    l_t(0,count1,count2) = 1
		    ! Calculate |l_0|^2
		    l_tsq(0,count1) = l_tsq(0,count1) + 1
		  else
		    l_t(0,count1,count2) = 0
		  endif
		endif
	      end do !m2
	    end do !s2
	  end do !m1
	end do !s1

100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.EQ.-1) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	call calc_com

	! Construct the neighbour list at the current time
	count1 = 0
	do s1=sp1,sp1end
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1
	    l_tsq(nframes,count1) = 0
	    c1x=comx(s1,m1)
	    c1y=comy(s1,m1)
	    c1z=comz(s1,m1)
	    count2 = 0
	    do s2=sp2,sp2end
	      do m2=1,s_nmols(s2)
		count2 = count2 + 1
		if ((s1.EQ.s2).AND.(m1.EQ.m2)) then
		  l_t(nframes,count1,count2) = 0
		else
		  call pbc(comx(s2,m2),comy(s2,m2),comz(s2,m2),c1x,c1y,c1z,tx,ty,tz)
		  rij=sqrt( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )
		  ! Heaviside function
		  if (rij.LT.cagecut) then
		    l_t(nframes,count1,count2) = 1
		    ! Calculate |l_t|^2 for origin summing
		    l_tsq(nframes,count1) = l_tsq(nframes,count1) + 1
		  else
		    l_t(nframes,count1,count2) = 0
		  endif
		endif
	      end do !m2
	    end do !s2
	  end do !m1
	end do !s1

	! Update cage correlation function origins
	do n=1,nframes
	  ! Evaluate : n_i_out = |l_0|^2 - l_0*l_t
	  !   Our origin (*_0) is defined by the frame loop above (n)
	  !   |l_0|^2 calculated in previous loop when constructing l_t(t,,) == l_t(n,)
	  !   l_0*l_t is calculated now: l_0t = l_t(n,) * l_t(nframes,)  [our origin * current frame]
	  count1 = 0
	  do s1=sp1,sp1end
	    do m1=1,s_nmols(s1)
	      count1 = count1 + 1
	      count2 = 0
	      l_0t = 0
	      do s2=sp2,sp2end
		 do m2=1,s_nmols(s2)
		  count2 = count2 + 1
		  l_0t = l_0t + l_t(n-1,count1,count2) * l_t(nframes,count1,count2)
		 end do !m2
	      end do !s2
	      ! Sum data for this origin.
	      n_i_out = l_tsq(n-1,count1) - l_0t
	      ! Heaviside function again:
	      if (n_i_out.LE.c) then
	        ccor(nframes - (n-1)) = ccor(nframes - (n-1)) + 1.0
	      end if
	    end do !m1
	  end do !s1
	  ! Increase origin normalisation counter
	  ccorn(nframes - (n-1)) = ccorn(nframes - (n-1)) + 1
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

	! Take averages of the data.
	! Write out time and natural logarithm of ccor (plus populations)
	open(unit=9,file="cagecor.dat",form="formatted",status="replace")
	do n=1,nframes
	  ccor(n) = ccor(n) / ccorn(n)
	  if (sp1.EQ.0) then
	    count1 = totmols
	  else
	    count1 = s_nmols(sp1)
	  end if
	  if (sp1.EQ.sp2) then
	    ccor(n) = ccor(n) / (count1-1)
	  else
	    ccor(n) = ccor(n) / count1
	  end if
	  write(9,"(3F10.6)") tstep*n, log(ccor(n)), ccor(n)
	end do
	close(9)

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program cagecor

