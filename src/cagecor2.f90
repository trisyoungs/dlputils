!	** cagecor2 **
!	Calculate the cage correlation function for species in the system
!	(Rabani, Gezelter, Berne - J. Chem. Phys, 1997, 107, 6867)

	program cagecor2
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,sp1,sp2,m1,m2,baselen,nframes,success,nargs,c, l_0t, n_i_out, maxnbrs
	integer :: framestodo, count, sp1end, sp2end, totmols, s1, s2, count1, count2, length
	integer :: ncentral, nouter, pos, tlast, t0, tn, nbr1, nbr2
	integer, allocatable :: nbrs(:,:,:), ccorn(:)
	integer :: iargc
	real*8, allocatable :: ccor(:)
	real*8 :: c1x,c1y,c1z,tx,ty,tz,rij,cagecut,cagecutsq,deltat

	nargs = iargc()
	if (nargs.LT.10) then
	  write(0,"(A)") "Usage : cagecor2 <HISTORYfile> <OUTPUTfile> <centresp> <outersp> <cagecutoff> <c> <maxnbrs> <length> <framestodo> <deltat>"
	  write(0,"(A)") "	 Set <centresp> or <outersp> to zero to specify 'any molecule'"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(I4)") sp1
	call getarg(4,temp); read(temp,"(I4)") sp2
	call getarg(5,temp); read(temp,"(F10.8)") cagecut
	cagecutsq = cagecut*cagecut
	call getarg(6,temp); read(temp,"(I6)") c
	call getarg(7,temp); read(temp,"(I6)") maxnbrs
	call getarg(8,temp); read(temp,"(I6)") length
	call getarg(9,temp); read(temp,"(I6)") framestodo
	call getarg(10,temp); read(temp,"(F10.8)") deltat

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
	! Now, read in the history header so that we have cell()
	if (readheader().EQ.-1) goto 799

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="ccresults."
	  baselen=10
	else
	  basename=hisfile(1:baselen)
	endif

	totmols = 0
	do n=1,nspecies
	  totmols = totmols + s_nmols(n)
	end do
	if (sp1.EQ.0) then
	  sp1 = 1
	  sp1end = nspecies
	  ncentral = totmols
	else
	  sp1end = sp1
	  ncentral = s_nmols(sp1)
	endif
	if (sp2.EQ.0) then
	  sp2 = 1
	  sp2end = nspecies
	  nouter = totmols
	else
	  sp2end = sp2
	  nouter = s_nmols(sp2)
	endif

	! Allocate arrays
	allocate (nbrs(length,ncentral,0:maxnbrs))
	allocate (ccor(0:length-1))
	allocate (ccorn(0:length-1))

	nbrs = 0
	ccor = 0.0
	ccorn = 0
	
100	nframes=0
	pos = 0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	! Calculate molecule centres of mass positions
	call calc_com

	! Determine 'position' in array and clear current data there
	pos = mod(nframes,length)+1
	nbrs(pos,:,:) = 0

	! Construct the neighbour list at the current time
	count1 = 0
	do s1=sp1,sp1end
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1
	    c1x=comx(s1,m1)
	    c1y=comy(s1,m1)
	    c1z=comz(s1,m1)
	    count2 = 0
	    do s2=sp2,sp2end
	      do m2=1,s_nmols(s2)
		count2 = count2 + 1
		if ((s1.EQ.s2).AND.(m1.EQ.m2)) cycle
		call pbc(comx(s2,m2),comy(s2,m2),comz(s2,m2),c1x,c1y,c1z,tx,ty,tz)
		tx = tx - c1x
		ty = ty - c1y
		tz = tz - c1z
		rij = tx*tx + ty*ty + tz*tz
		! Heaviside function
		if (rij.lt.cagecutsq) then
		  nbrs(pos,count1,0) = nbrs(pos,count1,0) + 1
		  nbrs(pos,count1,nbrs(pos,count1,0)) = count2
		endif
	      end do !m2
	    end do !s2
	  end do !m1
	end do !s1

	! Update cage correlation function origins
	if (nframes.lt.length) then
	  t0 = 1
	  tlast = nframes
	else
	  t0 = pos+1
	  if (t0.gt.length) t0 = t0-length
	  tlast = length
	end if

	do tn=1,tlast
	  n = tn-t0
	  if (n.lt.0) n = n + length
	  ! Evaluate : n_i_out = |l_0|^2 - l_0*nbrs
	  !   Our origin (*_0) is defined by the frame loop above (n)
	  !   |l_0|^2 calculated in previous loop when constructing nbrs(t,,) == nbrs(n,)
	  !   l_0*nbrs is calculated now: l_0t = nbrs(n,) * nbrs(nframes,)  [our origin * current frame]
	  count1 = 0
	  do s1=sp1,sp1end
	    do m1=1,s_nmols(s1)
	      count1 = count1 + 1
	      ! Compare lists of molecules between time origin t0 and current position 'n'
	      l_0t = 0
	      do nbr1=1,nbrs(t0,count1,0)
		do nbr2=1,nbrs(tn,count1,0)
		  if (nbrs(t0,count1,nbr1).eq.nbrs(tn,count1,nbr2)) l_0t = l_0t + 1
		end do
	      end do
	      ! Sum data for this origin.
	      n_i_out = nbrs(t0,count1,0) - l_0t
	      ! Heaviside function again:
	      if (n_i_out.le.c) ccor(n) = ccor(n) + 1.0
	    end do !m1
	  end do !s1
	  ! Increase origin normalisation counter
	  ccorn(n) = ccorn(n) + 1
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
	resfile=basename(1:baselen)//"cage"//CHAR(48+sp1)//CHAR(48+sp2)
	open(unit=9,file=resfile,form="formatted",status="replace")
	do n=0,length-1
	  ccor(n) = ccor(n) / ccorn(n) / ncentral
	  write(9,"(3F10.6,i10)") deltat*n, log(ccor(n)), ccor(n), ccorn(n)
	end do
	close(9)

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program cagecor2

