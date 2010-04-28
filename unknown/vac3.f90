!	** vac3 **
!	Calculate the velocity autocorrelation tensor over a restricted set of molecules

	program vacf3
	use dlprw; use utility
	implicit none
	character*80 :: hisfile1, hisfile2, dlpoutfile,basename,resfile,headerfile
	character*20 :: temp
	integer :: i,j,n1,n2,m,nframes,success,nargs,baselen,sp,spatom, vec(3)
	integer :: count1, selector, newflag, bin, finish, length, igeom(5)
	integer :: iargc
	integer, allocatable :: aflag(:,:)
	real*8, allocatable :: vx(:,:), vy(:,:), vz(:,:), vac(:), acc(:), vacpart(:,:), comvx(:,:), comvy(:,:), comvz(:,:)
	real*8, allocatable :: point(:)
	real*8 :: tx,ty,tz,scalar,deltat, geom(5), totmass, dist, vec1(3), vec2(3), centroid(2)
	real*8 :: xx, yy, zz, xy, xz, yz, mag, tempr
	logical :: testrun = .false., reverselogic = .false., altheader = .false., restart = .false.

	nargs = iargc()
	if (nargs.lt.8) then
	  write(0,"(a)") "Usage : vac3 <HISfile 1> <HISfile 2> <DLP OUTPUTfile> <headerfile [0 for none]> <delta t> <length> <species> <atom> <[-]selector> [extra data]"
	  stop
	end if
	call getarg(1,hisfile1)
	call getarg(2,hisfile2)
	call getarg(3,dlpoutfile)
	call getarg(4,headerfile)
	if (headerfile.ne."0") then
	  altheader = .true.
	  write(0,*) "Alternative header file supplied."
	end if
	call getarg(5,temp); read(temp,"(F10.6)") deltat
	call getarg(6,temp); read(temp,"(I6)") length
	if (length.lt.0) then
	  testrun = .true.
	  write(0,*) "Test run - selected comv's will be output."
	  open(unit=20,file='selectedpoints.xyz',status='replace',form='formatted')
	  length = abs(length)
	end if
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
	j = 10		! Position of first argument
	do i=j,nargs
	  temp = "                    "
	  call getarg(i,temp)
	  do n1=1,20
	    if (temp(n1:n1).eq.".") exit
	  end do
	  if (n1.lt.20) read(temp,"(f10.6)") geom(i-j+1)
	  if (n1.gt.20) read(temp,"(i10)") igeom(i-j+1)
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
	  case (6); write(0,"(a,'(',i2,'/',i2,'->',i2,')')") "ZOrient - Positive vectors of atoms ",igeom(1),igeom(2),igeom(3)
	  	    write(0,"(a,f9.4,a,f9.4)") "             between z = ",geom(4)," and ",geom(5)
		    vec(1) = igeom(1)
		    vec(2) = igeom(2)
		    vec(3) = igeom(3)
	  case (7); write(0,"(a,'(',i2,'/',i2,'->',i2,')')") "ZOrient - Negative vectors of atoms ",igeom(1),igeom(2),igeom(3)
	  	    write(0,"(a,f9.4,a,f9.4)") "             between z = ",geom(4)," and ",geom(5)
		    vec(1) = igeom(1)
		    vec(2) = igeom(2)
		    vec(3) = igeom(3)
	  case (8); write(0,"(a,f9.4,a,f9.4)") "ZHover - All molecules between z = ",geom(1)," and ",geom(2)
                    write(0,"(a,f9.4,a,f9.4,f9.4)") "             and within xy = ",geom(3)," of point ",geom(4), geom(5)
                    geom(3) = geom(3) * geom(3)
	end select

	if (outinfo(dlpoutfile,1).EQ.-1) goto 798

	allocate (aflag(length*2,natms))
	allocate (vac(0:length*2))
	allocate (vacpart(6,0:length*2))
	allocate (acc(0:length*2))
	allocate (vx(length*2,natms))
	allocate (vy(length*2,natms))
	allocate (vz(length*2,natms))
	allocate (comvx(length*2,s_nmols(sp)))
	allocate (comvy(length*2,s_nmols(sp)))
	allocate (comvz(length*2,s_nmols(sp)))
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
	
	! Read in data from restart file (if there is one)
	open(unit=20,file='vacfacc',form="unformatted",status="old",err=55)
	write(0,*) "Reading data from VACF accumulation file..."
	read(20,err=55) vac
	read(20,err=55) vacpart
	read(20,err=55) acc
	read(20) aflag
	restart = .true.
	close(20)
	goto 60

55	write(0,*) "No previous restart file 'vacfacc' found (or read error). Starting from zero..."
	vac = 0.0
	vacpart = 0.0
	acc = 0.0
	aflag = 0
	restart = .false.
	close(20)

	! To reduce memory overhead only store 'length' velocities in memory at any one time
	! Open first file and read 'length' frames from it
60	write(0,*) "Reading history file 1...."
	if (altheader) then
          call openhis(headerfile,10)
          if (readheader().EQ.-1) goto 799
          close(dlpun_his)
          call openhis(hisfile1,10)
	else 
	  call openhis(hisfile1,10)
	  if (readheader().EQ.-1) goto 799
	end if

100	nframes=0
!	Once we reach nframes=length, swap over to the other history file
101	if (nframes.eq.length) then
	  close(dlpun_his)
	  if (altheader) then
            call openhis(headerfile,10)
            if (readheader().EQ.-1) goto 799
            close(dlpun_his)
            call openhis(hisfile2,10)
	  else 
	    call openhis(hisfile2,10)
	    if (readheader().EQ.-1) goto 799
	  end if
	end if

	success=readframe()
	if (success.ne.0) goto 799  ! End of file encountered, or file error....

	nframes = nframes+1
	if (mod(nframes,100).EQ.0) write(0,"(i)") nframes

	! Store all velocities
	vx(nframes,:) = xvel(:)
	vy(nframes,:) = yvel(:)
	vz(nframes,:) = zvel(:)

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
	    ! Should also be normalised here, but we only require the sign of the z component
	    point(n1) = -(vec1(3) + vec2(3))
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
	      n1 = igeom(2)
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
	      n1 = igeom(2)
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
	      if ((point(j).gt.0.0).and.(zz.gt.geom(4)).and.(zz.lt.geom(5))) newflag = 1
	    case (7)	! Negative vectors within z-limits
	      if ((point(j).lt.0.0).and.(zz.gt.geom(4)).and.(zz.lt.geom(5))) newflag = 1
	    case (8)	! ZHover: within geom(1) < z < geom(2), and xy=geom(3) of point geom(4),geom(5)
	      if ((zz.gt.geom(1)).and.(zz.lt.geom(2))) then
		tx = xx - geom(4)
		ty = yy - geom(5)
		dist = tx*tx + ty*ty
		if (dist.lt.geom(3)) newflag = 1
	      end if
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

	  !If doing a restart, we only set aflags for the second 'lenght' frames
	  ! Set aflag() taking account of the previous state of the atom
	  if ((.not.restart).or.(restart.and.(nframes.gt.length))) then
	    if (nframes.eq.1) then
	      aflag(1,j) = newflag
	    else if (newflag.eq.1) then
	      ! Molecule is in region. Where was it before now?
	      if (aflag(nframes-1,j).le.1) then
	        aflag(nframes,j) = 1
	      else
	        aflag(nframes,j) = 2
	      end if
	    else
	      ! Not (or no longer) in region
	      if (aflag(nframes-1,j).eq.0) then
	        aflag(nframes,j) = 0
	      else
	        aflag(nframes,j) = 2
	      end if
	    end if
	  end if

	  ! For testrun, write out selected points
	  if (testrun) then
	    if (aflag(nframes,j).eq.1) write(20,"(a,3f12.5)") "H",xx,yy,zz
	  end if
	end do

	! Calculate COMV if necessary
	if (spatom.lt.1) then
	  ! Clear comv* for this position
	  comvx(nframes,:) = 0.0
	  comvy(nframes,:) = 0.0
	  comvz(nframes,:) = 0.0
	  i = s_start(sp)
	  do m=1,s_nmols(sp)
	    totmass = 0.0
	    do j=i,i+s_natoms(sp)-1
	      totmass = totmass + mass(j)
	      comvx(nframes,m) = comvx(nframes,m) + mass(j)*vx(nframes,j)
	      comvy(nframes,m) = comvy(nframes,m) + mass(j)*vy(nframes,j)
	      comvz(nframes,m) = comvz(nframes,m) + mass(j)*vz(nframes,j)
	    end do
	    !if (m.eq.1) write(0,*) "Frame ",nframes, "comvx(mol1) = ",comvx(nframes,1)
	    comvx(nframes,m) = comvx(nframes,m) / totmass
	    comvy(nframes,m) = comvy(nframes,m) / totmass
	    comvz(nframes,m) = comvz(nframes,m) / totmass
	    i = i + s_natoms(sp)
	  end do
	endif

	if (nframes.eq.length*2) goto 301
	goto 101

	! Accumulate VACF.
301	write(0,*) "Accumulating vacf..."

	! TEST - Allow molecules to re-enter regions and accumulate again...
	do i=1,nframes
	  do n1=1,s_nmols(sp)
	    if (aflag(i,n1).eq.2) aflag(i,n1) = 1
	  end do
	end do
	
	! Loop over target species/atoms and calculate only the consecutive part of the VACF where the
	! molecules were in the selected region for the first time.
	! Do up to 'length' velocities in succession, starting from arbitrary point 'pos' in the array
	do n1=1,length
	  i = 1
	  if (spatom.gt.0) i = s_start(sp) - 1 + spatom
	  if (mod(n1,100).eq.0) write(0,*) "Calculating VACF - current position = ",n1

	  do m=1,s_nmols(sp)
	  
	    ! Only accumulate if this molecule is selected in the origin frame n1
	    if (aflag(n1,i).ne.1) then
	      if (spatom.lt.1) then
	        i = i + 1
	      else
	        i = i + s_natoms(sp)
	      end if
	      cycle
	    end if

	    ! Consider all points in first 'length' frames with themselves and the second 'length' frames

	    do n2=n1,n1+length

	!write(0,*) "Considering n1,n2",n1,n2
	      ! Distance limit between frames is 'length'
	      if (n2-n1.ge.length) exit

	      ! If aflag() <> 1, we're done with this atom / molecule
	      if (aflag(n2,i).ne.1) exit

	      ! If aflag() == 0 the molecule is not (yet) in the region of interest
	      !if (aflag(n2,i).eq.0) cycle

	      ! Calculate VACF to time t0(n1) -> t1(n2)
	      ! Calculate products
	      if (spatom.lt.1) then
	        xx = comvx(n1,i)*comvx(n2,i)
	        yy = comvy(n1,i)*comvy(n2,i)
	        zz = comvz(n1,i)*comvz(n2,i)
	        xy = comvx(n1,i)*comvy(n2,i) + comvy(n1,i)*comvx(n2,i)
	        xz = comvx(n1,i)*comvz(n2,i) + comvz(n1,i)*comvx(n2,i)
	        yz = comvy(n1,i)*comvz(n2,i) + comvz(n1,i)*comvy(n2,i)
	      else
	        xx = vx(n1,i)*vx(n2,i)
	        yy = vy(n1,i)*vy(n2,i)
	        zz = vz(n1,i)*vz(n2,i)
	        xy = vx(n1,i)*vy(n2,i) + vy(n1,i)*vx(n2,i)
	        xz = vx(n1,i)*vz(n2,i) + vz(n1,i)*vx(n2,i)
	        yz = vy(n1,i)*vz(n2,i) + vz(n1,i)*vy(n2,i)
	      end if

	      ! Calculate full VACF
	      scalar = xx + yy + zz

	      ! Calculate bin and accumulate
	      bin = n2 - n1
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
	end do   ! vacf accumulation

	goto 801

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended before 'length*2' was fulfilled..."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",length*2,")"
	goto 801
800	write(0,*) "Framestodo was fulfilled."
801	write(0,*) ""

	! Ascertain length of basename....
	baselen=-1
	do i=80,1,-1
	  if (hisfile1(i:i).eq.".") then
	    baselen=i
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="vacresults."
	  baselen=11
	else
	  basename=hisfile1(1:baselen)
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
	  write(20,"(8f14.6, e20.6)") n1*deltat, vac(n1) / acc(n1), (vacpart(m,n1) / acc(n1),m=1,6), acc(n1)
	end do
	close(20)

	! Write restart data
	open(unit=20,file='vacfacc',form="unformatted",status="replace")
	write(20) vac
	write(20) vacpart
	write(20) acc
	write(20) aflag
	close(20)

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program vacf3
