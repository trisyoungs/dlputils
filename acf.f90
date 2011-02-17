!	** acf **
!	Calculate the specified autocorrelation tensor, optionally over a restricted set of molecules

	program acfunc
	use dlprw; use utility
	implicit none
	integer, parameter :: VELOCITY=1, MOLDIPOLE=2, SYSDIPOLE=3
	character*80 :: hisfile1, dlpoutfile,basename,resfile,headerfile
	character*20 :: temp
	integer :: i,j,n1,n2,m,nframes,success,nargs,baselen,sp,spatom, pos, lastpos
	integer :: count1, bin, n1finish, mfinish, length, acftype = 0, framestodo = -1
	integer :: iargc
	real*8, allocatable :: acf(:), accum(:), acfpart(:,:), qx(:,:), qy(:,:), qz(:,:)
	real*8 :: tx,ty,tz,scalar,deltat, totmass, dist, vec(3)
	real*8 :: xx, yy, zz, xy, xz, yz
	logical :: altheader = .false.

	nargs = iargc()
	if (nargs.ne.8) then
	  write(0,"(a)") "Usage : acf <HISfile> <DLP OUTPUTfile> <headerfile [0 for none]> <functype=velocity,moldipole,sysdipole> <delta t> <length> <species> <nframes>"
	  stop
	end if
	call getarg(1,hisfile1)
	call getarg(2,dlpoutfile)
	call getarg(3,headerfile)
	if (headerfile.ne."0") then
	  altheader = .true.
	  write(0,*) "Alternative header file supplied."
	end if
	call getarg(4,temp)
	if (temp.eq."velocity") acftype = VELOCITY
	if (temp.eq."moldipole") acftype = MOLDIPOLE
	if (temp.eq."sysdipole") acftype = SYSDIPOLE
	if (acftype.eq.0) stop "Invalid correlation function requested - options are 'velocity' or 'dipole'."
	call getarg(5,temp); read(temp,"(F10.6)") deltat
	call getarg(6,temp); read(temp,"(I6)") length
	call getarg(7,temp); read(temp,"(I6)") sp
	call getarg(8,temp); read(temp,"(I6)") framestodo
	
	if (acftype.eq.VELOCITY) write(0,*) "Calculating velocity autocorrelation function."
	if (acftype.eq.MOLDIPOLE) write(0,*) "Calculating molecular dipole autocorrelation function."
	if (acftype.eq.SYSDIPOLE) write(0,*) "Calculating system dipole autocorrelation function."
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798

	allocate (acf(0:length))
	allocate (acfpart(6,0:length))
	allocate (accum(0:length))
	if (acftype.eq.SYSDIPOLE) then
	  allocate (qx(length,1))
	  allocate (qy(length,1))
	  allocate (qz(length,1))
	else
	  allocate (qx(length,s_nmols(sp)))
	  allocate (qy(length,s_nmols(sp)))
	  allocate (qz(length,s_nmols(sp)))
	end if

	acf = 0.0
	acfpart = 0.0
	accum = 0.0
	qx = 0.0
	qy = 0.0
	qz = 0.0
	
	! To reduce memory overhead only store 'length' velocities in memory at any one time
60	write(0,*) "Reading history file...."
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
	pos=0
	lastpos=-1
101	success=readframe()
	if (success.ne.0) goto 799  ! End of file encountered, or file error....

	nframes = nframes+1
	if (mod(nframes,100).EQ.0) write(0,"(i)") nframes

	! Calculate array position and (wrapped) last position
	pos = mod(nframes-1,length) + 1
	lastpos = pos - 1
	if (lastpos.eq.0) lastpos = length

	! Calculate quantities for correlation function
	qx(pos,:) = 0.0
	qy(pos,:) = 0.0
	qz(pos,:) = 0.0
	if (acftype.eq.VELOCITY) then
	  call calc_com
	  i = s_start(sp)
	  do m=1,s_nmols(sp)
	    totmass = 0.0
	    do j=i,i+s_natoms(sp)-1
	      totmass = totmass + mass(j)
	      qx(pos,m) = qx(pos,m) + mass(j)*xvel(j)
	      qy(pos,m) = qy(pos,m) + mass(j)*yvel(j)
	      qz(pos,m) = qz(pos,m) + mass(j)*zvel(j)
	    end do
	    !if (m.eq.1) write(0,*) "Frame ",pos, "comvx(mol1) = ",comvx(pos,1)
	    qx(pos,m) = qx(pos,m) / totmass
	    qy(pos,m) = qy(pos,m) / totmass
	    qz(pos,m) = qz(pos,m) / totmass
	    i = i + s_natoms(sp)
	  end do
	else if (acftype.eq.MOLDIPOLE) then
	  ! Take all atom positions relative to first atom in molecule
	  i = s_start(sp)
	  do m=1,s_nmols(sp)
	    do j=i,i+s_natoms(sp)-1
	      ! Get mim position of this atom with first
	      call pbc(xpos(j),ypos(j),zpos(j),xpos(i),ypos(i),zpos(i),vec(1),vec(2),vec(3))
	      qx(pos,m) = qx(pos,m) + charge(j)*vec(1)
	      qy(pos,m) = qy(pos,m) + charge(j)*vec(2)
	      qz(pos,m) = qz(pos,m) + charge(j)*vec(3)
	    end do
	    i = i + s_natoms(sp)
	  end do
	else if (acftype.eq.SYSDIPOLE) then
	  ! Take all atom positions relative to first atom in molecule
	  i = s_start(sp)
	  do m=1,s_nmols(sp)
	    do j=i,i+s_natoms(sp)-1
	      ! Get mim position of this atom with first
	      call pbc(xpos(j),ypos(j),zpos(j),xpos(i),ypos(i),zpos(i),vec(1),vec(2),vec(3))
	      qx(pos,1) = qx(pos,1) + charge(j)*vec(1)
	      qy(pos,1) = qy(pos,1) + charge(j)*vec(2)
	      qz(pos,1) = qz(pos,1) + charge(j)*vec(3)
	    end do
	    i = i + s_natoms(sp)
	  end do
	endif

	! Accumulate VACF.
	! Loop over target species/atoms and calculate only the consecutive parts of the VACF
	! Do up to 'length' velocities in succession, starting from arbitrary point 'pos' in the array
	!write(0,*) "Calculating VACF - current frame = ",nframes, "pos = ",pos
	mfinish = s_nmols(sp)
	if (acftype.eq.SYSDIPOLE) mfinish = 1
	do m=1,mfinish

	  !write(0,"(i6,a,i6)") m,"/",s_nmols(sp)
	  ! 'pos' contains the position of the most recently stored velocities.
	  ! Loop over all other stored velocities, calculating scalars with this point.

	  ! Step through the velocities *backwards* from the last stored values
	  n1finish = max(pos-length+1,pos-nframes+1)
	  do n1=pos,n1finish,-1

	    ! 'Fold' n1 into the region 1-length
	    n2 = n1
	    if (n2.lt.1) n2 = n2 + length

	    ! Calculate ACF to time t0(n1) -> t1(n2)
	    ! Calculate products
	    xx = qx(n2,m)*qx(pos,m)
	    yy = qy(n2,m)*qy(pos,m)
	    zz = qz(n2,m)*qz(pos,m)
	    xy = qx(n2,m)*qy(pos,m) + qy(n2,m)*qx(pos,m)
	    xz = qx(n2,m)*qz(pos,m) + qz(n2,m)*qx(pos,m)
	    yz = qy(n2,m)*qz(pos,m) + qz(n2,m)*qy(pos,m)

	    ! Calculate full VACF
	    scalar = xx + yy + zz

	    ! Calculate bin - 'pos' is the last position read, so 'n1' *must* be less than pos when we calculate the bin
	    bin = pos - n2
	    if (bin.lt.0) bin = bin + length
	    !if (bin.eq.0) write(55,"(a,3i4,a,i3,a,4f12.4)") "Acc ",pos,n2,m, "bin ",bin, " scalar ",scalar,comvx(pos,i),comvy(pos,i),comvz(pos,i)
	    acf(bin) = acf(bin) + scalar
	    accum(bin) = accum(bin) + 1.0
	    ! Calculate partial VACFs
	    acfpart(1,bin) = acfpart(1,bin) + xx
	    acfpart(2,bin) = acfpart(2,bin) + yy
	    acfpart(3,bin) = acfpart(3,bin) + zz
	    acfpart(4,bin) = acfpart(4,bin) + xy
	    acfpart(5,bin) = acfpart(5,bin) + xz
	    acfpart(6,bin) = acfpart(6,bin) + yz

	  end do

	end do

	if (nframes.eq.framestodo) goto 801
	goto 101

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended before 'length*2' was fulfilled..."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",nframes,")"
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
	if (acftype.eq.VELOCITY) resfile=basename(1:baselen)//"vacf"//CHAR(48+sp)
	if (acftype.eq.MOLDIPOLE) resfile=basename(1:baselen)//"mdacf"//CHAR(48+sp)
	if (acftype.eq.SYSDIPOLE) resfile=basename(1:baselen)//"sdacf"//CHAR(48+sp)
	open(unit=20,file=resfile,form="formatted",status="replace")

	! Write the data...
	write(20,"(9a14)") "#t","acf","xx","yy","zz","xy","xz","yz","acc"
	do n1=0,length-1
	  write(20,"(8f14.6, e20.6)") n1*deltat, acf(n1) / accum(n1), (acfpart(m,n1) / accum(n1),m=1,6), accum(n1)
	end do
	close(20)

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program acfunc
