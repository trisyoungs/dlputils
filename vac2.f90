!	** vac2 **
!	Calculate the velocity autocorrelation tensor over a restricted set of molecules

	program vacf2
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile,headerfile
	character*20 :: temp
	integer :: i,j,n1,n1start,n1finish,m,nframes,success,nargs,baselen,sp,spatom, vec(3)
	integer :: framestodo, count1, selector, newflag, bin, finish, length, pos, lastpos, igeom(6), n2
	integer :: iargc
	integer, allocatable :: aflag(:,:)
	real*8, allocatable :: vx(:,:), vy(:,:), vz(:,:), vac(:), acc(:), vacpart(:,:), comvx(:,:), comvy(:,:), comvz(:,:)
	real*8, allocatable :: point(:)
	real*8 :: tx,ty,tz,scalar,deltat, geom(6), totmass, dist, vec1(3), vec2(3)
	real*8 :: xx, yy, zz, xy, xz, yz, mag, pvec(3)
	logical :: testrun = .false., reverselogic = .false., altheader = .false.

	nargs = iargc()
	if (nargs.lt.8) then
	  write(0,"(a)") "Usage : vac2 <HISTORYfile> <OUTPUTfile> <headerfile [0 for none]> <delta t> <framestodo> <length> <sp> <atom> <[-]selector> [extra data]"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,headerfile)
	if (headerfile.ne."0") then
	  altheader = .true.
	  write(0,*) "Alternative header file supplied."
	end if
	call getarg(4,temp); read(temp,"(F10.6)") deltat
	call getarg(5,temp); read(temp,"(I6)") framestodo
	if (framestodo.lt.0) then
	  testrun = .true.
	  write(0,*) "Test run - selected comv's will be output."
	  open(unit=20,file='selectedpoints.xyz',status='replace',form='formatted')
	  framestodo = abs(framestodo)
	end if
	call getarg(6,temp); read(temp,"(I6)") length
	call getarg(7,temp); read(temp,"(I6)") sp
	call getarg(8,temp); read(temp,"(I6)") spatom
	call getarg(9,temp); read(temp,"(I6)") selector
	if (selector.lt.0) then
	  selector = abs(selector)
	  write(0,*) "Selector logic will be *reversed*."
	  reverselogic = .true.
	end if
	if ((sp.lt.1).or.(spatom.lt.0)) then
	  write(0,*) "Negative number found for sp or spatom."
	  stop
	end if
	
	geom = 0.0
	do i=10,nargs
	  temp = "                    "
	  call getarg(i,temp)
	  do n1=1,20
	    if (temp(n1:n1).eq.".") exit
	  end do
	  if (n1.lt.20) read(temp,"(f10.6)") geom(i-9)
	  if (n1.gt.20) read(temp,"(i10)") igeom(i-9)
	end do
	if (nargs.gt.8) write(0,*) "Extra data is:", geom
	if (spatom.lt.1) write(0,*) "Using COM velocity instead of specific particle velocity."
	if ((selector.eq.6).and.(spatom.gt.0)) stop "Pointing vectors require the COM velocity to be used."

	! Print info on selector and do any geom() adjustments
	write(0,*) "Defined selector is:"
	select case (selector)
	  case (0); write(0,"(a)") "None - all molecules of species will be considered."
	  case (1); write(0,"(a,f9.4,a,f9.4)") "Zslice - All molecules between z = ",geom(1)," and ",geom(2)
	  case (2); write(0,"(a,f9.4,a,f9.4)") "Yslice - All molecules between y = ",geom(1)," and ",geom(2)
	  case (3); write(0,"(a,f9.4,a,f9.4)") "Xslice - All molecules between x = ",geom(1)," and ",geom(2)
	  case (4); write(0,"(a,f9.4,a,f2.1)") "Proximity - All molecules within ",geom(1)," of species ",geom(2)
		    geom(1) = geom(1) * geom(1)
	  case (5); write(0,"(a,f9.4,a,i2)") "ZProximity - All molecules within ",geom(1)," of species ",igeom(2)
	  	    write(0,"(a,f9.4,a,f9.4)") "             and between z = ",geom(3)," and ",geom(4)
		    geom(1) = geom(1) * geom(1)
	  case (6); write(0,"(a,'(',i2,'/',i2,'->',i2,')')") "ZOrient - Half positive vectors of atoms ",igeom(1),igeom(2),igeom(3)
	  	    write(0,"(a,f9.4,a,f9.4)") "             between z = ",geom(4)," and ",geom(5)
		    vec(1) = igeom(1)
		    vec(2) = igeom(2)
		    vec(3) = igeom(3)
	  case (7); write(0,"(a,'(',i2,'/',i2,'->',i2,')')") "ZOrient - Half negative vectors of atoms ",igeom(1),igeom(2),igeom(3)
	  	    write(0,"(a,f9.4,a,f9.4)") "             between z = ",geom(4)," and ",geom(5)
		    vec(1) = igeom(1)
		    vec(2) = igeom(2)
		    vec(3) = igeom(3)
	  case (8); write(0,"(a,f9.4,a,f9.4)") "XYZbox - All molecules between z = ",geom(1)," and ",geom(2)
	  	    write(0,"(a,f9.4,a,f9.4)") "XYZbox - Within distance of ",geom(3)," from X= ",geom(4)
	  	    write(0,"(a,f9.4,a,f9.4)") "XYZbox - Within distance of ",geom(5)," from Y= ",geom(6)
	end select

	! Open and check the files...
	if (altheader) then
          call openhis(headerfile,10)
          if (readheader().EQ.-1) goto 799
          close(dlpun_his)
          call openhis(hisfile,10)
	else 
	  call openhis(hisfile,10)
	  if (readheader().EQ.-1) goto 799
	end if

	if (outinfo(dlpoutfile,1).EQ.-1) goto 798

	allocate (aflag(length,natms))
	allocate (vac(0:length))
	allocate (vacpart(6,0:length))
	allocate (acc(0:length))
	allocate (vx(length,natms))
	allocate (vy(length,natms))
	allocate (vz(length,natms))
	allocate (comvx(length,s_nmols(sp)))
	allocate (comvy(length,s_nmols(sp)))
	allocate (comvz(length,s_nmols(sp)))
	allocate (point(s_nmols(sp)))

	vac = 0.0
	vacpart = 0.0
	acc = 0.0
	aflag = 0
	vx = 0.0
	vy = 0.0
	vz = 0.0
	comvx = 0.0
	comvy = 0.0
	comvz = 0.0
	
	! To reduce memory overhead only store 'length' velocities in memory at any one time
	! To prevent copying of memory to maintain linear array of velocities, use mod(nframes,length)
100	nframes=0
	pos=0
	lastpos = -1
101	success=readframe()
	if (success.ne.0) goto 799  ! End of file encountered, or file error....

	nframes = nframes+1
	if (mod(nframes,100).EQ.0) write(0,"(i6)") nframes

	! Calculate array position and (wrapped) last position
	pos = mod(nframes-1,length) + 1
	lastpos = pos - 1
	if (lastpos.eq.0) lastpos = length
	!write(0,*) pos,lastpos

	! Store all velocities
	vx(pos,:) = xvel(:)
	vy(pos,:) = yvel(:)
	vz(pos,:) = zvel(:)

	! Perform 'selection' on atoms / COMs (flag those that are / were in the correct region)
	! aflag values: 0 = not in region, never has been
	! 		1 = currently in region for the first time
	!		2 = (not) in region but has been previously so ignore
	! Calculate COM of molecules if necessary
	if ((spatom.lt.1).or.(selector.ge.4)) call calc_com

	! Calculate pointing vectors if necessary
	if ((selector.eq.6).or.(selector.eq.7)) then
	  i = s_start(sp)
	  do n1=1,s_nmols(sp)
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
	    point(n1) = -pvec(3)

	    i = i + s_natoms(sp)
	  end do
	end if

	if (spatom.lt.1) then
	  finish = s_nmols(sp)
	else
	  finish = natms
	end if

	do j=1,finish
	  ! Get point to consider
	  if (spatom.lt.1) then
	    xx = comx(sp,j)
	    yy = comy(sp,j)
	    zz = comz(sp,j)
	  else
	    xx = xpos(j)
	    yy = ypos(j)
	    zz = zpos(j)
	  end if
	  newflag = 0

	  ! Do selection of atoms 
	  select case (selector)
	    case (0)	! No geometric selection (do all atoms)
	      newflag = 1
	    case (1)	! Between slice in z-direction (geom(1) and geom(2)
	      if ((zz.gt.geom(1)).and.(zz.lt.geom(2))) newflag = 1
	    case (2)	! Between slice in y-direction (geom(1) and geom(2)
	      if ((yy.gt.geom(1)).and.(yy.lt.geom(2))) newflag = 1
	    case (3)	! Between slice in x-direction (geom(1) and geom(2)
	      if ((xx.gt.geom(1)).and.(xx.lt.geom(2))) newflag = 1
	    case (4)	! Within geom(1) angstroms of a molecule in species geom(2)
	      n1 = nint(geom(2))
	      do m=1,s_nmols(n1)
	        call pbc(xx,yy,zz,comx(n1,m),comy(n1,m),comz(n1,m),tx,ty,tz)
		tx = tx - comx(n1,m)
		ty = ty - comy(n1,m)
		tz = tz - comz(n1,m)
		dist = tx*tx + ty*ty + tz*tz
		if (dist.lt.geom(1)) newflag = 1
		if (newflag.eq.1) exit
	      end do
	    case (5)	! Within geom(1) angstroms of a molecule in species geom(2) and z-limits
	      n1 = nint(geom(2))
	      do m=1,s_nmols(n1)
	        call pbc(xx,yy,zz,comx(n1,m),comy(n1,m),comz(n1,m),tx,ty,tz)
		tx = tx - comx(n1,m)
		ty = ty - comy(n1,m)
		tz = tz - comz(n1,m)
		dist = tx*tx + ty*ty + tz*tz
	!write(0,"(i4,8f8.3)") m,xx,yy,zz,comx(n1,m),comy(n1,m), comz(n1,m),dist,geom(1)
		if ((dist.lt.geom(1)).and.(zz.gt.geom(3)).and.(zz.lt.geom(4))) newflag = 1
		if (newflag.eq.1) exit
	      end do
	    case (6)	! Positive vectors within z-limits
	      if ((point(j).gt.0.5).and.(zz.gt.geom(4)).and.(zz.lt.geom(5))) newflag = 1
	    case (7)	! Negative vectors within z-limits
	      if ((point(j).lt.-0.5).and.(zz.gt.geom(4)).and.(zz.lt.geom(5))) newflag = 1
	    case (8)	! Boxed in specified volume
	      call pbc(xx,yy,zz,geom(4),geom(6),0.5*(geom(1)+geom(2)),tx,ty,tz)
	      tx = dabs(tx - geom(4))
	      ty = dabs(ty - geom(6))
	      newflag = 1
	!write(0,*) xx,yy,zz
	      if ((tz.lt.geom(1)).or.(tz.gt.geom(2))) newflag = 0
	      if (tx.gt.geom(3)) newflag = 0
	      if (ty.gt.geom(5)) newflag = 0
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
	    aflag(pos,j) = 1
	  else
	    ! Not in region
	    aflag(pos,j) = 0
	  end if

	  ! For testrun, write out selected points
	  if (testrun) then
	    if (aflag(pos,j).eq.1) write(20,"(a,3f12.5)") "H",xx,yy,zz
	  end if
	end do

	! Calculate COMV if necessary
	if (spatom.lt.1) then
	  ! Clear comv* for this position
	  comvx(pos,:) = 0.0
	  comvy(pos,:) = 0.0
	  comvz(pos,:) = 0.0
	  i = s_start(sp)
	  do m=1,s_nmols(sp)
	    totmass = 0.0
	    do j=i,i+s_natoms(sp)-1
	      totmass = totmass + mass(j)
	      comvx(pos,m) = comvx(pos,m) + mass(j)*vx(pos,j)
	      comvy(pos,m) = comvy(pos,m) + mass(j)*vy(pos,j)
	      comvz(pos,m) = comvz(pos,m) + mass(j)*vz(pos,j)
	    end do
	    !if (m.eq.1) write(0,*) "Frame ",nframes, "comvx(mol1) = ",comvx(pos,1)
	    comvx(pos,m) = comvx(pos,m) / totmass
	    comvy(pos,m) = comvy(pos,m) / totmass
	    comvz(pos,m) = comvz(pos,m) / totmass
	    i = i + s_natoms(sp)
	  end do
	endif

	! Accumulate VACF.
	! Loop over target species/atoms and calculate only the consecutive parts of the VACF
	! Do up to 'length' velocities in succession, starting from arbitrary point 'pos' in the array
	i = 1
	if (spatom.gt.0) i = s_start(sp) - 1 + spatom
	!write(0,*) "Calculating VACF - current frame = ",nframes, "pos = ",pos
	do m=1,s_nmols(sp)

	  !write(0,"(i6,a,i6)") m,"/",s_nmols(sp)
	  ! 'pos' contains the position of the most recently stored velocities.
	  ! Loop over all other stored velocities, calculating scalars with this point.

	  ! Only accumulate if the reference position is inside the region (aflag(pos,i) = 1)
	  if (aflag(pos,i).eq.0) then
	    if (spatom.lt.1) then
	      i = i + 1
	    else
	      i = i + s_natoms(sp)
	    end if
	    cycle
	  endif

	  ! Step through the velocities *backwards* from the last stored values
	  n1finish = max(pos-length+1,pos-nframes+1)
	  do n1=pos,n1finish,-1

	    ! 'Fold' n1 into the region 1-length
	    n2 = n1
	    if (n2.lt.1) n2 = n2 + length
	    ! If aflag() == 0 the molecule is not in the region of interest (anymore) so go to next molecule
	    if (aflag(n2,i).eq.0) exit

	    ! Calculate VACF to time t0(n1) -> t1(n2)
	    ! Calculate products
	    if (spatom.lt.1) then
	      xx = comvx(n2,i)*comvx(pos,i)
	      yy = comvy(n2,i)*comvy(pos,i)
	      zz = comvz(n2,i)*comvz(pos,i)
	      xy = comvx(n2,i)*comvy(pos,i) + comvy(n2,i)*comvx(pos,i)
	      xz = comvx(n2,i)*comvz(pos,i) + comvz(n2,i)*comvx(pos,i)
	      yz = comvy(n2,i)*comvz(pos,i) + comvz(n2,i)*comvy(pos,i)
	    else
	      xx = vx(n2,i)*vx(pos,i)
	      yy = vy(n2,i)*vy(pos,i)
	      zz = vz(n2,i)*vz(pos,i)
	      xy = vx(n2,i)*vy(pos,i) + vy(n2,i)*vx(pos,i)
	      xz = vx(n2,i)*vz(pos,i) + vz(n2,i)*vx(pos,i)
	      yz = vy(n2,i)*vz(pos,i) + vz(n2,i)*vy(pos,i)
	    end if

	    ! Calculate full VACF
	    scalar = xx + yy + zz

	    ! Calculate bin - 'pos' is the last position read, so 'n1' *must* be less than pos when we calculate the bin
	    bin = pos - n2
	    if (bin.lt.0) bin = bin + length
	    !if (bin.eq.0) write(55,"(a,3i4,a,i3,a,4f12.4)") "Acc ",pos,n2,m, "bin ",bin, " scalar ",scalar,comvx(pos,i),comvy(pos,i),comvz(pos,i)
	    vac(bin) = vac(bin) + scalar
	    acc(bin) = acc(bin) + 1.0
	    ! Calculate partial VACFs
	    vacpart(1,bin) = vacpart(1,bin) + xx
	    vacpart(2,bin) = vacpart(2,bin) + yy
	    vacpart(3,bin) = vacpart(3,bin) + zz
	    vacpart(4,bin) = vacpart(4,bin) + xy
	    vacpart(5,bin) = vacpart(5,bin) + xz
	    vacpart(6,bin) = vacpart(6,bin) + yz

	  end do

	  if (spatom.lt.1) then
	    i = i + 1
	  else
	    i = i + s_natoms(sp)
	  end if

	end do

	! Next frame?
	if (nframes.eq.framestodo) goto 801
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
	do i=80,1,-1
	  if (hisfile(i:i).eq.".") then
	    baselen=i
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="vacresults."
	  baselen=11
	else
	  basename=hisfile(1:baselen)
	endif

	! Open output file
	if (spatom.lt.1) then
	  resfile=basename(1:baselen)//"vacfcom"//CHAR(48+sp)
	else
	  resfile=basename(1:baselen)//"vacf"//CHAR(48+sp)//"a"//CHAR(48+spatom)
	end if
	open(unit=20,file=resfile,form="formatted",status="replace")

	! Write the data...
	write(20,"(9a14)") "#t","vacf","xx","yy","zz","xy","xz","yz","acc"
	do n1=0,length-1
	  write(20,"(8f14.6, f20.6)") n1*deltat, vac(n1) / acc(n1), (vacpart(m,n1) / acc(n1),m=1,6), acc(n1)
	end do

	close(20)

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program vacf2

