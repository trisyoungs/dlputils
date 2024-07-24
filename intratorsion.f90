!	** torsion **
!	Calculate torsion angle histogram

	program torsionhist
	use dlprw; use utility
	implicit none
	real*8, parameter :: pi = 3.14159265358979d0, radcon = 57.29577951d0
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*2 :: a1, a2, a3, a4
	character*20 :: temp
	integer :: status,  nargs, success, baselen
	integer :: aoff1,n,m,s1,s2,bin,nframes,numadded,sp,atom1,atom2,atom3,atom4
	integer :: iargc, ntorsionbins
	real*8 :: i(3), j(3), k(3), l(3), v(3), vecji(3), vecjk(3), veckj(3), veckl(3), xp1(3), xp2(3)
	real*8 :: norm, dp, torsion, rij, tx, ty, tz, ktx, kty, ktz, mag1, mag2, torsionbin
	real*8, allocatable :: ijkl(:)
	real*8 :: dist, magnitude, dotproduct

	torsionbin = 1.0
	ntorsionbins = 180.0 / torsionbin

	nargs = iargc()
	if (nargs.ne.7) then
	  write(0,*) "Usage : intratorsion <HISTORYfile> <OUTPUTfile> <sp> <i> <j> <k> <l>"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(i4)") sp
	call getarg(4,temp); read(temp,"(i4)") atom1
	call getarg(5,temp); read(temp,"(i4)") atom2
	call getarg(6,temp); read(temp,"(i4)") atom3
	call getarg(7,temp); read(temp,"(i4)") atom4

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	write(0,"(A,I5)") "Target species is ",sp
	write(0,"(a,4(i3))") "Calculating torsion angle from atoms ",atom1,atom2,atom3,atom4
	allocate(ijkl(-ntorsionbins:ntorsionbins),stat=status); if (status.GT.0) stop "Allocation error for ijkl()"
	
	! Initialise the arrays...
	ijkl = 0.0

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

	  ! Grab coordinates of the atoms
	  i(1) = xpos(aoff1+atom1-1)
	  i(2) = ypos(aoff1+atom1-1)
	  i(3) = zpos(aoff1+atom1-1)

	  j(1) = xpos(aoff1+atom2-1)
	  j(2) = ypos(aoff1+atom2-1)
	  j(3) = zpos(aoff1+atom2-1)

	  k(1) = xpos(aoff1+atom3-1)
	  k(2) = ypos(aoff1+atom3-1)
	  k(3) = zpos(aoff1+atom3-1)

	  l(1) = xpos(aoff1+atom4-1)
	  l(2) = ypos(aoff1+atom4-1)
	  l(3) = zpos(aoff1+atom4-1)

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
	  ! Calculate sign
	  dp = dotproduct(xp1, veckl)
	    if (dp.lt.0) then
	      torsion = -(torsion+1.0)
	      bin = int(torsion / torsionbin)
	      if (bin.lt.-ntorsionbins) bin = -ntorsionbins
	    else
	      bin = int(torsion / torsionbin)
	      if (bin.ge.ntorsionbins) bin = ntorsionbins-1
	    end if
	  ijkl(bin) = ijkl(bin) + 1

	  ! Global counter
	  numadded = numadded+1

	  aoff1 = aoff1 + s_natoms(sp)
	end do   ! End main loop over molecules

	if (nframes.EQ.1) write(0,*) "numadded:=",numadded
	if (nframes.EQ.1) then
	  write(0,"(A,I2,A,I4,A)") "Histogram calculated for ",sp," over ",s_nmols(sp)," molecules."
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
	ijkl = ijkl / nframes

	a1 = char(48+atom1/10)//char(48+MOD(atom1,10))
	a2 = char(48+atom2/10)//char(48+MOD(atom2,10))
	a3 = char(48+atom3/10)//char(48+MOD(atom3,10))
	a4 = char(48+atom4/10)//char(48+MOD(atom4,10))

	! Write torsion histogram ijkl
	resfile=basename(1:baselen)//a1//"-"//a2//"-"//a3//"-"//a4//".tors"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=-ntorsionbins,ntorsionbins-1
	  write(9,"(f9.3,2x,f12.4)") (n+0.5)*torsionbin,ijkl(n)
	end do
	close(9)

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end program torsionhist

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

