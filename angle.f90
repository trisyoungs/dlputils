!	** angle **
!	Calculate angle between vectors on different molecules

	program anglecalc
	use dlprw; use utility
	implicit none
	real*8, parameter :: pi = 3.14159265358979d0
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*2 :: a1, a2
	character*20 :: temp
	integer :: status,  nargs, success, baselen
	integer :: nanglebins,ndistbins,aoff1,aoff2,n,m,s1,s2,bin,bin2,nframes,numadded,sp,atom1,atom2
	integer :: iargc
	real*8 :: i1(3), i2(3), j1(3), j2(3), ij1(3), ij2(3), v(3),dp, distbin, anglebin, norm, angle, dist
	real*8, allocatable :: ijj(:), ijj_jj(:,:), iji(:), jij(:)

	anglebin=0.5	! In degrees
	distbin=0.1	! In Angstroms

	nargs = iargc()
	if (nargs.ne.5) then
	  write(0,*) "Usage : angle <HISTORYfile> <OUTPUTfile> <sp> <atom1> <atom2>"
	  write(0,*) "             --- Geometry is i1-j1...j2-i2"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(i4)") sp
	call getarg(4,temp); read(temp,"(i4)") atom1
	call getarg(5,temp); read(temp,"(i4)") atom2

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	nanglebins = 180.0 / anglebin
	ndistbins = cell(1) / distbin
	write(0,"(A,F6.3,A,F6.3,A)") "Using binwidths of ",anglebin,"and",distbin," Degrees / Angstroms (Angle / Distance)"
	write(0,"(A,I5,A)") "There will be ",nanglebins," angle bins."
	write(0,"(A,I5,A)") "There will be ",ndistbins," distance bins."
	write(0,"(A,I5)") "Target species is ",sp
	write(0,"(a,i2,a,i2)") "Calculating angle between vectors from atoms ",atom1," and ",atom2
	
	allocate(ijj(nanglebins),stat=status); if (status.GT.0) stop "Allocation error for ijj()"
	allocate(iji(nanglebins),stat=status); if (status.GT.0) stop "Allocation error for iji()"
	allocate(jij(nanglebins),stat=status); if (status.GT.0) stop "Allocation error for jij()"
	allocate(ijj_jj(nanglebins,ndistbins),stat=status); if (status.GT.0) stop "Allocation error for ijj_jj()"

	! Initialise the arrays...
	ijj = 0.0
	iji = 0.0
	jij = 0.0
	ijj_jj = 0.0

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 120  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	
	aoff1 = s_start(sp)
	do s1 = 1,s_nmols(sp)     ! Loop over all molecules in species

	  ! Grab coordinates of the first vector (i1->j1)
	  i1(1) = xpos(aoff1+atom1-1)
	  i1(2) = ypos(aoff1+atom1-1)
	  i1(3) = zpos(aoff1+atom1-1)

	  j1(1) = xpos(aoff1+atom2-1)
	  j1(2) = ypos(aoff1+atom2-1)
	  j1(3) = zpos(aoff1+atom2-1)

	  ! PBC atoms
	  call pbc(j1(1),j1(2),j1(3),i1(1),i1(2),i1(3),ij1(1),ij1(2),ij1(3))
	  ! Store PBC position of j
	  j1 = ij1
	  ! Calculate vector and normalise
	  ij1 = ij1 - i1
	  dist = sqrt(ij1(1)*ij1(1) + ij1(2)*ij1(2) + ij1(3)*ij1(3))
	  ij1 = ij1 / dist

	  aoff2 = s_start(sp)
	  do s2 = 1,s_nmols(sp)

	    if (s1.eq.s2) then
	      aoff2 = aoff2 + s_natoms(sp)
	      cycle
	    end if

	    ! Grab coordinates of second vector (i2->j2)
	    i2(1) = xpos(aoff2+atom1-1)
	    i2(2) = ypos(aoff2+atom1-1)
	    i2(3) = zpos(aoff2+atom1-1)

	    j2(1) = xpos(aoff2+atom2-1)
	    j2(2) = ypos(aoff2+atom2-1)
	    j2(3) = zpos(aoff2+atom2-1)

	    ! PBC i2 w.r.t. i1, and j2 w.r.t. i2
	    call pbc(i2(1),i2(2),i2(3),i1(1),i1(2),i1(3),ij2(1),ij2(2),ij2(3))
	    i2 = ij2
	    call pbc(j2(1),j2(2),j2(3),i2(1),i2(2),i2(3),ij2(1),ij2(2),ij2(3))
	    j2 = ij2
	    ! Calculate vector i2->j2 and normalise
	    ij2 = ij2 - i2
	    dist = sqrt(ij2(1)*ij2(1) + ij2(2)*ij2(2) + ij2(3)*ij2(3))
	    ij2 = ij2 / dist

	    !
	    ! Analysis
	    !
	    ! 1) Angle i1-j1-j2
	    ! Get vector j1->j2 and normalise
	    v = j2 - j1
	    dist = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
	    v = v / dist
	    dp = ij1(1)*v(1) + ij1(2)*v(2) + ij1(3)*v(3)
	    angle = safeAngle(dp)
	    bin = int(angle * (1.0 / anglebin)) + 1
	    ijj(bin) = ijj(bin)+1
	    ! write(0,*) dp,angle,dist

	    ! 2) Distance j1-j2 (for matrix histogram)
	    ! 'dist' calculated above already has the distance, so can just use this...
	    bin2 = int(dist * (1.0/distbin)) + 1
	    if (bin2.le.ndistbins) ijj_jj(bin,bin2) = ijj_jj(bin,bin2) + 1

	    ! 3) Angle i1-j1-i2
	    ! Get vector j1->i2 and normalise
	    v = i2 - j1
	    dist = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
	    v = v / dist
	    dp = ij1(1)*v(1) + ij1(2)*v(2) + ij1(3)*v(3)
	    angle = safeAngle(dp)
	    bin = int(angle * (1.0 / anglebin)) + 1
	    iji(bin) = iji(bin)+1

	    ! 4) Angle j1-i1/2-j2
	    dp = ij1(1)*ij2(1) + ij1(2)*ij2(2) + ij1(3)*ij2(3)
	    angle = safeAngle(dp)
	    bin = int(angle * (1.0 / anglebin)) + 1
	    jij(bin) = jij(bin)+1

	    ! Global counter
	    numadded = numadded+1

	    aoff2 = aoff2 + s_natoms(sp)
	  end do
	  aoff1 = aoff1 + s_natoms(sp)
	end do   ! End main loop over all atoms of species1.

	if (nframes.EQ.1) write(0,*) "numadded:=",numadded
	if (nframes.EQ.1) then
	  write(0,"(A,I2,A,I4,A)") "PRDF of atoms about ",sp," : averaged over ",s_nmols(sp)," molecules."
	end if

	! Next frame
	goto 101
	! Work on the results now to get the proper RDF
120	write(0,*) "Finished."
	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended prematurely!"
	goto 801
800	write(0,*) "End of unformatted HISTORY file found."
801	write(0,*) ""

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).EQ.".") THEN
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) THEN
	  basename="angleresults."
	  baselen=11
	else
	  basename=hisfile(1:baselen)
	endif

	! Normalise arrays.
	! Use expected total, not numadded, for uniformity (since some points will never be binned)
	norm = nframes * s_nmols(sp) * (s_nmols(sp) - 1)
	ijj = ijj / nframes
	iji = iji / nframes
	jij = jij / nframes
	ijj_jj = ijj_jj / nframes

	a1 = char(48+atom1/10)//char(48+MOD(atom1,10))
	a2 = char(48+atom2/10)//char(48+MOD(atom2,10))

	! Write histogram ijj
	resfile=basename(1:baselen)//a1//"-"//a2//".ijj"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,nanglebins
	  write(9,"(F10.3,3x,e12.6)") anglebin*(n-0.5),ijj(n)
	end do
	close(9)

	! Write histogram iji
	resfile=basename(1:baselen)//a1//"-"//a2//".iji"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,nanglebins
	  write(9,"(F10.3,3x,e12.6)") anglebin*(n-0.5),iji(n)
	end do
	close(9)

	! Write histogram jij
	resfile=basename(1:baselen)//a1//"-"//a2//".jij"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,nanglebins
	  write(9,"(F10.3,3x,e12.6)") anglebin*(n-0.5),jij(n)
	end do
	close(9)

	! Write matrix histogram
	resfile=basename(1:baselen)//a1//"-"//a2//".ijj_jj"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,nanglebins
	  do m=1,ndistbins
	    write(9,"(f10.3,2x,f6.3,2x,e12.6)") anglebin*(n-0.5),distbin*(m-0.5),ijj_jj(n,m)
	  end do
	  write(9,*) ""
	end do
	close(9)

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end program anglecalc

