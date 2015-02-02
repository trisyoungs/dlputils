!	** msd2 **
!	Calculate the mean square displacement over a subset of molecules (COM)
!	Reads all frames at once, and calculates msd per molecule (slow)

	program msd2
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,m,m1,nframes,success,nargs,baselen,i,j
	integer :: framestodo, count1, sp, selector, ncaught
	integer :: iargc, igeom(6), newflag, vec(3)
	integer, allocatable :: aflag(:,:)
	logical :: reverselogic = .false.
	real*8, allocatable :: msd(:,:), point(:), finalmsd(:,:)
	real*8, allocatable :: x(:,:), y(:,:), z(:,:)
	integer, allocatable :: orn(:)
	real*8 :: tx,ty,tz,rij2,deltat, geom(6), pvec(3), vec1(3), vec2(3), mag, xx, yy, zz, dist, dx, dy, dz

	nargs = iargc()
	if (nargs.lt.6) then
	  write(0,"(A)") "Usage : msd <HISTORYfile> <OUTPUTfile> <delta t> <framestodo> <sp> <[-]selector> [extra data]"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(F10.6)") deltat
	call getarg(4,temp); read(temp,"(I6)") framestodo
	call getarg(5,temp); read(temp,"(I6)") sp
	call getarg(6,temp); read(temp,"(I6)") selector
	if (selector.lt.0) then
	  selector = abs(selector)
	  write(0,*) "Selector logic will be *reversed*."
	  reverselogic = .true.
	end if

	geom = 0.0
	do i=7,nargs
	  temp = "                    "
	  call getarg(i,temp)
	  do n=1,20
	    if (temp(n:n).eq.".") exit
	  end do
	  if (n.lt.20) read(temp,"(f10.6)") geom(i-6)
	  if (n.gt.20) read(temp,"(i10)") igeom(i-6)
	end do
	if (nargs.gt.6) write(0,*) "Extra data is:", geom

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
	    goto 50
	  endif
	end do
50     if (baselen.EQ.-1) then
	  basename="rdfresults."
	  baselen=11
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=11,file=basename(1:baselen)//"caught",form="formatted",status="replace")
	open(unit=9,file=basename(1:baselen)//"msd",form="formatted",status="replace")
	! Print info on selector and do any geom() adjustments
	write(0,*) "Defined selector is:"
	select case (selector)
	  case (0); write(9,"(a)") "# None - all molecules of species will be considered."
	  case (1); write(9,"(a,f9.4,a,f9.4)") "# Zslice - All molecules between z = ",geom(1)," and ",geom(2)
	  case (2); write(9,"(a,f9.4,a,f9.4)") "# Yslice - All molecules between y = ",geom(1)," and ",geom(2)
	  case (3); write(9,"(a,f9.4,a,f9.4)") "# Xslice - All molecules between x = ",geom(1)," and ",geom(2)
	  case (4); write(9,"(a,f9.4,a,f2.1)") "# Proximity - All molecules within ",geom(1)," of species ",geom(2)
		    geom(1) = geom(1) * geom(1)
	  case (5); write(9,"(a,f9.4,a,i2)") "# ZProximity - All molecules within ",geom(1)," of species ",igeom(2)
	  	    write(9,"(a,f9.4,a,f9.4)") "#              and between z = ",geom(3)," and ",geom(4)
		    geom(1) = geom(1) * geom(1)
	  case (6); write(9,"(a,'(',i2,'/',i2,'->',i2,')')") "# ZOrient - Positive half vectors of atoms ",igeom(1),igeom(2),igeom(3)
	  	    write(9,"(a,f9.4,a,f9.4)") "#              between z = ",geom(4)," and ",geom(5)
		    vec(1) = igeom(1)
		    vec(2) = igeom(2)
		    vec(3) = igeom(3)
	  case (7); write(9,"(a,'(',i2,'/',i2,'->',i2,')')") "# ZOrient - Negative half vectors of atoms ",igeom(1),igeom(2),igeom(3)
	  	    write(9,"(a,f9.4,a,f9.4)") "#              between z = ",geom(4)," and ",geom(5)
		    vec(1) = igeom(1)
		    vec(2) = igeom(2)
		    vec(3) = igeom(3)
	end select

	allocate (aflag(framestodo,s_nmols(sp)))
	allocate (point(s_nmols(sp)))

	allocate (orn(0:framestodo))
	allocate (msd(0:framestodo,4))
	allocate (finalmsd(0:framestodo,4))
	allocate (x(framestodo,s_nmols(sp)))
	allocate (y(framestodo,s_nmols(sp)))
	allocate (z(framestodo,s_nmols(sp)))

	orn = 0
	aflag = 0
	msd = 0.0
	finalmsd = 0.0
	
	! Read in all COM coordinates and determine selections
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	! Store folded COM coordinates
	call calc_com
	do m1=1,s_nmols(sp)
	  call dlpfold(comx(sp,m1),comy(sp,m1),comz(sp,m1),x(nframes,m1),y(nframes,m1),z(nframes,m1))
	  !x(nframes,m1) = comx(sp,m1)
	  !y(nframes,m1) = comy(sp,m1)
	  !z(nframes,m1) = comz(sp,m1)
	  !x(nframes,m1) = xpos(s_start(sp)+(m1-1)*s_natoms(sp))
	  !y(nframes,m1) = ypos(s_start(sp)+(m1-1)*s_natoms(sp))
	  !z(nframes,m1) = zpos(s_start(sp)+(m1-1)*s_natoms(sp))
	end do !m1

	! Perform 'selection' on COMs (flag those that are / were in the correct region)
	! aflag values: 0 = not in region
	! 		1 = currently in region

	ncaught = 0
	! Calculate pointing vectors if necessary
	if ((selector.eq.6).or.(selector.eq.7)) then
	  i = s_start(sp)
	  do n=1,s_nmols(sp)
	    ! Get mim unit vector position1 from first atom with third atom (in vec1)
	    call pbc(xpos(i+vec(1)-1),ypos(i+vec(1)-1),zpos(i+vec(1)-1), &
		& xpos(i+vec(3)-1),ypos(i+vec(3)-1),zpos(i+vec(3)-1), vec1(1), vec1(2), vec1(3))
	    vec1(1) = vec1(1) - xpos(i+vec(3)-1)
	    vec1(2) = vec1(2) - ypos(i+vec(3)-1)
	    vec1(3) = vec1(3) - zpos(i+vec(3)-1)
	    mag = sqrt(vec1(1)*vec1(1) + vec1(2)*vec1(2) + vec1(3)*vec1(3))
	    vec1(1) = vec1(1) / mag
	    vec1(2) = vec1(2) / mag
	    vec1(3) = vec1(3) / mag
	    
	    ! Get MIM unit vector position from second atom with third atom (in vec2)
	    call pbc(xpos(i+vec(2)-1),ypos(i+vec(2)-1),zpos(i+vec(2)-1), &
		& xpos(i+vec(3)-1),ypos(i+vec(3)-1),zpos(i+vec(3)-1), vec2(1), vec2(2), vec2(3))
	    vec2(1) = vec2(1) - xpos(i+vec(3)-1)
	    vec2(2) = vec2(2) - ypos(i+vec(3)-1)
	    vec2(3) = vec2(3) - zpos(i+vec(3)-1)
	    mag = sqrt(vec2(1)*vec2(1) + vec2(2)*vec2(2) + vec2(3)*vec2(3))
	    vec2(1) = vec2(1) / mag
	    vec2(2) = vec2(2) / mag
	    vec2(3) = vec2(3) / mag

	    ! Average of vectors vec1 and vec2 give the point which our pointing vector runs through (from O)
	    ! So, pointing vector is |t+p| which is the (H1-H2)->O vector
	    pvec(1) = vec1(1) + vec2(1)
	    pvec(2) = vec1(2) + vec2(2)
	    pvec(3) = vec1(3) + vec2(3)
	    mag = sqrt(pvec(1)*pvec(1) + pvec(2)*pvec(2) + pvec(3)*pvec(3))
	    pvec(1) = pvec(1) / mag
	    pvec(2) = pvec(2) / mag
	    pvec(3) = pvec(3) / mag

	    ! Should also be normalised here, but we only require the sign of the z component
	    point(n) = -pvec(3)
	    i = i + s_natoms(sp)
	  end do
	end if

	do j=1,s_nmols(sp)
	  ! Get point to consider
	  xx = x(nframes,j)
	  yy = y(nframes,j)
	  zz = z(nframes,j)
	  newflag = 0

	  ! Do selection of atoms 
	  select case (selector)
	    case (0)	! No geometric selection (do all atoms)
	      newflag = 1
	    case (1)	! Between slice in z-direction (geom(1) and geom(2)
	      if ((zz.gt.geom(1)).and.(zz.lt.geom(2))) newflag = 1
		!if (newflag.eq.0) write(0,*) "Molecule ",j," rejected as its folded z-COM is ",zz
	    case (2)	! Between slice in y-direction (geom(1) and geom(2)
	      if ((yy.gt.geom(1)).and.(yy.lt.geom(2))) newflag = 1
	    case (3)	! Between slice in x-direction (geom(1) and geom(2)
	      if ((xx.gt.geom(1)).and.(xx.lt.geom(2))) newflag = 1
	    case (4)	! Within geom(1) angstroms of a molecule in species geom(2)
	      n = nint(geom(2))
	      do m=1,s_nmols(n)
	        call pbc(xx,yy,zz,comx(n,m),comy(n,m),comz(n,m),tx,ty,tz)
		tx = tx - comx(n,m)
		ty = ty - comy(n,m)
		tz = tz - comz(n,m)
		dist = tx*tx + ty*ty + tz*tz
		if (dist.lt.geom(1)) newflag = 1
		if (newflag.eq.1) exit
	      end do
	    case (5)	! Within geom(1) angstroms of a molecule in species geom(2) and z-limits
	      n = nint(geom(2))
	      do m=1,s_nmols(n)
	        call pbc(xx,yy,zz,comx(n,m),comy(n,m),comz(n,m),tx,ty,tz)
		tx = tx - comx(n,m)
		ty = ty - comy(n,m)
		tz = tz - comz(n,m)
		dist = tx*tx + ty*ty + tz*tz
	!write(0,"(i4,8f8.3)") m,xx,yy,zz,comx(n1,m),comy(n1,m), comz(n1,m),dist,geom(1)
		if ((dist.lt.geom(1)).and.(zz.gt.geom(3)).and.(zz.lt.geom(4))) newflag = 1
		if (newflag.eq.1) exit
	      end do
	    case (6)	! Positive vectors within z-limits
	      if ((point(j).gt.0.5).and.(zz.gt.geom(4)).and.(zz.lt.geom(5))) newflag = 1
	    case (7)	! Negative vectors within z-limits
	      if ((point(j).lt.-0.5).and.(zz.gt.geom(4)).and.(zz.lt.geom(5))) newflag = 1
	    case default
	      write(0,*) "No criteria coded for selector ",selector
	      stop
	  end select
	  ! Use reverse selector logic if required
	  if (reverselogic) then
	    if (newflag.eq.1) then
	      newflag = 0
	    else
	      newflag = 1
	    end if
	  end if

	  ! Set aflag()
	  if (nframes.eq.1) then
	    aflag(1,j) = newflag
	  else if (newflag.eq.1) then
	    ! Molecule is in region.
	    aflag(nframes,j) = 1
	  else
	    ! Not in region
	    aflag(nframes,j) = 0
	  end if

	  ! Increase ncaught
	  if (newflag.eq.1) ncaught = ncaught + 1

	end do

	write(11,*) nframes, ncaught

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

	write(0,*) "Calculating MSD...."

	! Loop over molecules in the target species, accumulating the position delta as we go
	do m1=1,s_nmols(sp)

	  if (mod(m1,10).eq.0) write(0,*) m1

	  ! Loop over frames (outer)
	  do i=1,nframes

	    ! Skip if this molecule isn't selected
	    if (aflag(i,m1).eq.0) cycle

	    dx = 0.0
	    dy = 0.0
	    dz = 0.0

	    ! Loop over frames (inner)
	    do j=i+1,nframes

	      ! End this run of accumulation if the molecule leaves the region of interest
	      if (aflag(j,m1).ne.1) exit

	      ! Calculate position change between this frame and the last
	      call pbc(x(j,m1),y(j,m1),z(j,m1),x(j-1,m1),y(j-1,m1),z(j-1,m1),tx,ty,tz)
	      tx = tx - x(j-1,m1)
	      ty = ty - y(j-1,m1)
	      tz = tz - z(j-1,m1)
	      !if (tx.gt.5.0) write(0,"(a,i4,a,i5,a,f)") "X-delta for molecule ",m1," at frame ",j," seems rather large...",tx
	      !if (ty.gt.5.0) write(0,"(a,i4,a,i5,a,f)") "Y-delta for molecule ",m1," at frame ",j," seems rather large...",ty
	      !if (tz.gt.5.0) write(0,"(a,i4,a,i5,a,f)") "Z-delta for molecule ",m1," at frame ",j," seems rather large...",tz
	      dx = dx + tx
	      dy = dy + ty
	      dz = dz + tz

	      ! Calculate MSD
	      xx = dx*dx
	      yy = dy*dy
	      zz = dz*dz
	      rij2 = xx + yy + zz

	      ! Increase total msd_xyz
	      msd(j-i,1) = msd(j-i,1) + rij2
	      msd(j-i,2) = msd(j-i,2) + xx
	      msd(j-i,3) = msd(j-i,3) + yy
	      msd(j-i,4) = msd(j-i,4) + zz
	      orn(j-i) = orn(j-i) + 1

	    end do	!j

	  end do	!i


	end do	!m1

	
	! Average over number of molecules and origins
	do n=0,nframes

	  if (orn(n).ne.0) then
	    finalmsd(n,1) = msd(n,1) / real(orn(n))
	    finalmsd(n,2) = msd(n,2) / real(orn(n))
	    finalmsd(n,3) = msd(n,3) / real(orn(n))
	    finalmsd(n,4) = msd(n,4) / real(orn(n))
	  end if
	
	  ! Write the data...
	  write(9,"(8e15.8)") n*deltat, (finalmsd(n,m),m=1,4), real(orn(n))
	end do
	close(9)

	write(0,*) "Finished."
999	close(10)
	close(11)
	close(13)
	end program msd2

