!	** torsion **
!	Calculate torsion angle formed between vectors on molecules of a species
! 	Angles and distances are output also

	program torsioncalc
	use dlprw; use utility
	implicit none
	real*8, parameter :: pi = 3.14159265358979d0, radcon = 57.29577951d0
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*2 :: a1, a2, a3, a4
	character*20 :: temp
	integer :: status,  nargs, success, baselen
	integer :: ntorsionbins,nanglebins,ndistbins,aoff1,aoff2,n,m,s1,s2,bin,bin2,nframes,numadded,sp1,sp2,atom_i,atom_j,atom_k,atom_l
	integer :: iargc
	real*8 :: i(3), j(3), k(3), l(3), v(3), vecji(3), vecjk(3), veckj(3), veckl(3), xp1(3), xp2(3)
	real*8 :: dp, distbin, anglebin, torsionbin, torsion, angle, rij, tx, ty, tz, ktx, kty, ktz, mag1, mag2, cutoff
	real*8, allocatable :: ijk(:), jkl(:), ijkl(:), ijkl_jk(:,:), ijkl_ijk(:,:)
	real*8 :: dist, magnitude, dotproduct

	torsionbin=0.5	! In degrees
	anglebin=0.5	! In degrees
	distbin=0.1	! In Angstroms
	cutoff=-1.0

	nargs = iargc()
	if (nargs.lt.8) then
	  write(0,*) "Usage : intertorsion <HISTORYfile> <OUTPUTfile> <sp1> <i> <j> <sp2> <k> <l> [maxjk]"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(i4)") sp1
	call getarg(4,temp); read(temp,"(i4)") atom_i
	call getarg(5,temp); read(temp,"(i4)") atom_j
	call getarg(6,temp); read(temp,"(i4)") sp2
	call getarg(7,temp); read(temp,"(i4)") atom_k
	call getarg(8,temp); read(temp,"(i4)") atom_l
	if (nargs.eq.9) then
	  call getarg(9,temp)
	  read(temp,"(f8.4)") cutoff
	endif

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	ntorsionbins = 180.0 / torsionbin
	nanglebins = 180.0 / anglebin
	ndistbins = cell(1) / distbin
	write(0,"(A,F6.3,F6.3,F6.3)") "Distance, angle, torsion binwidths are ",distbin,anglebin,torsionbin
	write(0,"(A,I5,A)") "There will be ",nanglebins," angle bins."
	write(0,"(A,I5,A)") "There will be ",ndistbins," distance bins."
	if (cutoff.gt.0) write(0,"(A,F7.3)") "Enforcing j-k cutoff of ",cutoff
	write(0,"(A,i5,i5)") "Target species are ",sp1, sp2
	write(0,"(a,i2,a,i2,a,i2,a,i2,a)") "Calculating i(sp1,",atom_i,")-j(sp1,",atom_j,")-k(sp2,",atom_k,")-l(sp2,",atom_l,")"
	
	allocate(ijk(nanglebins),stat=status); if (status.GT.0) stop "Allocation error for ijk()"
	allocate(jkl(nanglebins),stat=status); if (status.GT.0) stop "Allocation error for jkl()"
	allocate(ijkl(-ntorsionbins:ntorsionbins),stat=status); if (status.GT.0) stop "Allocation error for ijkl()"
	allocate(ijkl_jk(-ntorsionbins:ntorsionbins,ndistbins),stat=status); if (status.GT.0) stop "Allocation error for ijkl_jk()"
	allocate(ijkl_ijk(-ntorsionbins:ntorsionbins,nanglebins),stat=status); if (status.GT.0) stop "Allocation error for ijkl_ijk()"

	! Initialise the arrays...
	ijk = 0.0
	jkl = 0.0
	ijkl = 0.0
	ijkl_jk = 0.0
	ijkl_ijk = 0.0

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
	numadded = 0
101	success=readframe()
	if (success.EQ.1) goto 120  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	
	aoff1 = s_start(sp1)
	do s1 = 1,s_nmols(sp1)     ! Loop over all molecules in species

	  ! Grab coordinates of the first vector (i->j)
	  i(1) = xpos(aoff1+atom_i-1)
	  i(2) = ypos(aoff1+atom_i-1)
	  i(3) = zpos(aoff1+atom_i-1)

	  j(1) = xpos(aoff1+atom_j-1)
	  j(2) = ypos(aoff1+atom_j-1)
	  j(3) = zpos(aoff1+atom_j-1)

	  aoff2 = s_start(sp2)
	  do s2 = 1,s_nmols(sp2)

	    if ((sp1.eq.sp2).and.(s1.eq.s2)) then
	      aoff2 = aoff2 + s_natoms(sp2)
	      cycle
	    end if

	    ! Grab coordinates of second vector (k->l)
	    k(1) = xpos(aoff2+atom_k-1)
	    k(2) = ypos(aoff2+atom_k-1)
	    k(3) = zpos(aoff2+atom_k-1)

	    l(1) = xpos(aoff2+atom_l-1)
	    l(2) = ypos(aoff2+atom_l-1)
	    l(3) = zpos(aoff2+atom_l-1)

	    !
	    ! Calculate vectors
	    !
	    ! Angle i-j-k
	    ! Vector j->i
	    call pbc(i(1),i(2),i(3),j(1),j(2),j(3),tx,ty,tz)
	    call getvector(tx,ty,tz,j(1),j(2),j(3),vecji(1),vecji(2),vecji(3))
	    ! Vector j->k
	    call pbc(k(1),k(2),k(3),j(1),j(2),j(3),ktx,kty,ktz)
	    call getvector(ktx,kty,ktz,j(1),j(2),j(3),vecjk(1),vecjk(2),vecjk(3))
	    ! Cutoff check
	    rij = magnitude(vecjk)
	    if ((cutoff.gt.0).and.(rij.gt.cutoff)) then
	      aoff2 = aoff2 + s_natoms(sp2)
	      cycle
	    end if
	    ! Angle j-k-l (mim w.r.t. k (mim j))
	    ! Vector k->j
	    veckj = -vecjk
	    ! Vector k->l
	    call pbc(l(1),l(2),l(3),ktx,kty,ktz,tx,ty,tz)
	    call getvector(tx,ty,tz,ktx,kty,ktz,veckl(1),veckl(2),veckl(3))

	    !
	    ! Calculate torsion angle
	    !
	    ! Calculate cross products and magnitudes
	    call crossproduct(vecji,vecjk,xp1)
	    mag1 = magnitude(xp1)
	    call crossproduct(veckj,veckl,xp2)
	    mag2 = magnitude(xp2)
	    ! Calculate dot product and angle...
	    dp = dotproduct(xp1, xp2) / (mag1 * mag2)
	    torsion = safeAngle(dp)
	    ! Calculate sign and bin
	    dp = dotproduct(xp1, veckl)
	    if (dp.lt.0) then
	      torsion = -(torsion+torsionbin)
	      bin = int(torsion / torsionbin)
	      if (bin.lt.-ntorsionbins) bin = -ntorsionbins
	    else
	      bin = int(torsion / torsionbin)
	      if (bin.ge.ntorsionbins) bin = ntorsionbins-1
	    end if
	    ijkl(bin) = ijkl(bin) + 1
	    ! write(0,*) dp,angle,dist

	    !
	    ! Calculate j1-j2 (j-k) distance to do matrix
	    !
	    rij = magnitude(vecjk)
	    bin2 = int(rij * (1.0/distbin)) + 1
	    if (bin2.le.ndistbins) ijkl_jk(bin,bin2) = ijkl_jk(bin,bin2) + 1
	    
	    !
	    ! Calculate angles i-j-k (i1-j1-j2) and j-k-l (j1-j2-i2)
	    !
	    ! i-j-k
	    dp = dotproduct(vecji,vecjk) / (magnitude(vecji) * magnitude(vecjk))
	    angle = safeAngle(dp)
	    bin2 = int(angle * (1.0/anglebin)) + 1
	    ijk(bin2) = ijk(bin2) + 1
	    ijkl_ijk(bin,bin2) = ijkl_ijk(bin,bin2) + 1
	    ! j-k-l
	    dp = dotproduct(veckj,veckl) / (magnitude(veckj) * magnitude(veckl))
	    angle = safeAngle(dp)
	    bin2 = int(angle * (1.0/anglebin)) + 1
	    jkl(bin2) = jkl(bin2) + 1


	    ! Global counter
	    numadded = numadded+1

	    aoff2 = aoff2 + s_natoms(sp2)
	  end do
	  aoff1 = aoff1 + s_natoms(sp1)
	end do   ! End main loop over all atoms of species1.

	if (nframes.EQ.1) write(0,*) "numadded:=",numadded

	! Next frame
	goto 101

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

	! Normalise arrays per frame
	ijk = ijk / nframes
	jkl = jkl / nframes
	ijkl = ijkl / nframes
	ijkl_jk = ijkl_jk / nframes
	ijkl_ijk = ijkl_ijk / nframes

	a1 = char(48+atom_i/10)//char(48+MOD(atom_i,10))
	a2 = char(48+atom_j/10)//char(48+MOD(atom_j,10))
	a3 = char(48+atom_k/10)//char(48+MOD(atom_k,10))
	a4 = char(48+atom_l/10)//char(48+MOD(atom_l,10))

	! Write distance histogram jk
	resfile=basename(1:baselen)//a1//"-"//a2//"-"//a3//"-"//a4//".jk"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,ndistbins
	  write(9,"(f7.3,2x,f12.8)") distbin*(n-0.5),sum(ijkl_jk(:,n))
	end do
	close(9)

	! Write angle histogram ijk
	resfile=basename(1:baselen)//a1//"-"//a2//"-"//a3//"-"//a4//".ijk"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,nanglebins
	  write(9,"(F7.3,2(3x,F12.8))") anglebin*(n-0.5), ijk(n), ijk(n)/sin(anglebin*(n-0.5)/radcon)
	end do
	close(9)

	! Write angle histogram jkl
	resfile=basename(1:baselen)//a1//"-"//a2//"-"//a3//"-"//a4//".jkl"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,nanglebins
	  write(9,"(F7.3,2(3x,F12.8))") anglebin*(n-0.5), jkl(n), jkl(n)/sin(anglebin*(n-0.5)/radcon)
	end do
	close(9)

	! Write torsion histogram ijkl
	resfile=basename(1:baselen)//a1//"-"//a2//"-"//a3//"-"//a4//".ijkl"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=-ntorsionbins,ntorsionbins-1
	  write(9,"(f10.3,2x,f12.8)") (n+0.5)*torsionbin,ijkl(n)
	end do
	close(9)

	! Write torsion/distance histogram ijkl_jk
	resfile=basename(1:baselen)//a1//"-"//a2//"-"//a3//"-"//a4//".ijkl_jk"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=-ntorsionbins,ntorsionbins-1
	  do m=1,ndistbins
	    write(9,"(f10.3,2x,f6.3,2x,f12.8)") (n+0.5)*torsionbin,distbin*(m-0.5),ijkl_jk(n,m)
	  end do
	  write(9,*) ""
	end do
	close(9)

	! Write torsion/angle histogram ijkl_ijk
	resfile=basename(1:baselen)//a1//"-"//a2//"-"//a3//"-"//a4//".ijkl_ijk"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=-ntorsionbins,ntorsionbins-1
	  do m=1,nanglebins
	    write(9,"(f10.3,2x,f6.3,2(3x,f12.8))") (n+0.5)*torsionbin,anglebin*(m-0.5), ijkl_ijk(n,m)/sin(anglebin*(n-0.5)/radcon)
	  end do
	  write(9,*) ""
	end do
	close(9)

	! Write torsion/angle histogram ijkl_ijk_norm
	resfile=basename(1:baselen)//a1//"-"//a2//"-"//a3//"-"//a4//".ijkl_ijk_norm"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=-ntorsionbins,ntorsionbins-1
	  do m=1,nanglebins
	    write(9,"(f10.3,2x,f10.3,2(3x,f12.8))") (n+0.5)*torsionbin,anglebin*(m-0.5), ijkl_ijk(n,m)
	  end do
	  write(9,*) ""
	end do
	close(9)

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end program torsioncalc

        real*8 function magnitude(vec)
	implicit none
        real*8 :: vec(3)
        magnitude=dsqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
        end function magnitude

        subroutine crossproduct(abc,xyz,result)
	implicit none
        real*8 :: abc(3), xyz(3), result(3)
        result(1)=abc(2)*xyz(3) - abc(3)*xyz(2)
        result(2)=abc(3)*xyz(1) - abc(1)*xyz(3)
        result(3)=abc(1)*xyz(2) - abc(2)*xyz(1)
        end subroutine crossproduct

        real*8 function dotproduct(abc,xyz)
	implicit none
        real*8 :: abc(3), xyz(3), result
        result=abc(1)*xyz(1) + abc(2)*xyz(2) + abc(3)*xyz(3)
        dotproduct = result
        end function dotproduct

