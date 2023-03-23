!	** dist2 **
!	Calculate the radial distribution function over a subset of molecules (COM)

	program dist2
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,m,m1,m2,nframes,success,nargs,baselen,i,j, bin, comatom = 0
	integer :: framestodo, count1, sp, selector
	integer :: iargc, igeom(6), newflag, vec(3), nbins
	integer, allocatable :: aflag(:)
	logical :: reverselogic = .false.
	real*8, allocatable :: hist(:), point(:), xx(:), yy(:), zz(:)
	real*8 :: tx,ty,tz, geom(6), pvec(3), vec1(3), vec2(3), mag, dist, binwidth

	nargs = iargc()
	if (nargs.lt.6) then
	  write(0,"(A)") "Usage : dist2 <HISTORYfile> <OUTPUTfile> <framestodo> <sp> <comatom or 0> <[-]selector> [extra data]"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(I6)") framestodo
	call getarg(4,temp); read(temp,"(I6)") sp
	call getarg(5,temp); read(temp,"(I6)") comatom
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

	open(unit=9,file=basename(1:baselen)//"dist",form="formatted",status="replace")
	! Print info on selector and do any geom() adjustments
	write(9,*) "# Defined selector was:"
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

	allocate (aflag(s_nmols(sp)))
	allocate (point(s_nmols(sp)))
	allocate (xx(s_nmols(sp)))
	allocate (yy(s_nmols(sp)))
	allocate (zz(s_nmols(sp)))

	binwidth = 0.025
	nbins = maxval(cell) / binwidth
	allocate (hist(0:nbins))

	aflag = 0
	hist = 0.0
	
	! Read in all COM coordinates and determine selections
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	! Store COM coordinates or comatom positions
	call calc_com
	if (comatom.lt.1) then
	  do m1=1,s_nmols(sp)
	    xx(m1) = comx(sp,m1)
	    yy(m1) = comy(sp,m1)
	    zz(m1) = comz(sp,m1)
	  end do !m1
	else
	  count1 = s_start(sp) - 1 + comatom
	  do m1=1,s_nmols(sp)
	    xx(m1) = xpos(count1)
	    yy(m1) = ypos(count1)
	    zz(m1) = zpos(count1)
	    count1 = count1 + s_natoms(sp)
	  end do !m1
	end if

	! Perform 'selection' on COMs (flag those that are / were in the correct region)
	! aflag values: 0 = not in region
	! 		1 = currently in region

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
	  newflag = 0

	  ! Do selection of atoms 
	  select case (selector)
	    case (0)	! No geometric selection (do all atoms)
	      newflag = 1
	    case (1)	! Between slice in z-direction (geom(1) and geom(2)
	      if ((zz(j).gt.geom(1)).and.(zz(j).lt.geom(2))) newflag = 1
		!if (newflag.eq.0) write(0,*) "Molecule ",j," rejected as its folded z-COM is ",zz
	    case (2)	! Between slice in y-direction (geom(1) and geom(2)
	      if ((yy(j).gt.geom(1)).and.(yy(j).lt.geom(2))) newflag = 1
	    case (3)	! Between slice in x-direction (geom(1) and geom(2)
	      if ((xx(j).gt.geom(1)).and.(xx(j).lt.geom(2))) newflag = 1
	    case (4)	! Within geom(1) angstroms of a molecule in species geom(2)
	      n = nint(geom(2))
	      do m=1,s_nmols(n)
	        call pbc(xx(j),yy(j),zz(j),comx(n,m),comy(n,m),comz(n,m),tx,ty,tz)
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
	        call pbc(xx(j),yy(j),zz(j),comx(n,m),comy(n,m),comz(n,m),tx,ty,tz)
		tx = tx - comx(n,m)
		ty = ty - comy(n,m)
		tz = tz - comz(n,m)
		dist = tx*tx + ty*ty + tz*tz
	!write(0,"(i4,8f8.3)") m,xx,yy,zz,comx(n1,m),comy(n1,m), comz(n1,m),dist,geom(1)
		if ((dist.lt.geom(1)).and.(zz(j).gt.geom(3)).and.(zz(j).lt.geom(4))) newflag = 1
		if (newflag.eq.1) exit
	      end do
	    case (6)	! Positive vectors within z-limits
	      if ((point(j).gt.0.5).and.(zz(j).gt.geom(4)).and.(zz(j).lt.geom(5))) newflag = 1
	    case (7)	! Negative vectors within z-limits
	      if ((point(j).lt.-0.5).and.(zz(j).gt.geom(4)).and.(zz(j).lt.geom(5))) newflag = 1
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
	  aflag(j) = newflag

	end do

	! Calculate histogram
	do m1=1,s_nmols(sp)
	  if (aflag(m1).eq.0) cycle
	  do m2=1,s_nmols(sp)
	    if (aflag(m2).eq.0) cycle
	    !call pbc(comx(sp,m2),comy(sp,m2),comz(sp,m2),comx(sp,m1),comy(sp,m1),comz(sp,m1),tx,ty,tz)
	    call pbc(xx(m2),yy(m2),zz(m2),xx(m1),yy(m1),zz(m1),tx,ty,tz)
	    tx = tx - xx(m1)
	    ty = ty - yy(m1)
	    tz = tz - zz(m1)
	    dist = sqrt(tx*tx + ty*ty + tz*tz)
	    bin = int(dist / binwidth)
	    if (bin.lt.nbins) hist(bin) = hist(bin) + 1
	  end do
	end do

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

	! Average over number of frames
	hist = hist / nframes
	
	do n=0,nbins
	  ! Write the data...
	  write(9,"(8e15.8)") n*binwidth, hist(n)
	end do
	close(9)

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program dist2

