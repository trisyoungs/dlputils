!	** intratorsion2 **
!	Calculate intramolecular torsion angle map

	program intratorsion2
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*2 :: a(8)
	character*20 :: temp
	integer :: status,  nargs, success, baselen
	integer :: aoff1,n,m,s1,s2,bin1,bin2,nframes,numadded,sp,t1(4),t2(4)
	integer :: iargc, ntorsionbins
	real*8 :: torsion, torsionbin
	real*8 :: psi1, psi2
	real*8, allocatable :: ijkl(:,:)

	torsionbin = 1.0
	ntorsionbins = 180 / torsionbin

	nargs = iargc()
	if (nargs.ne.11) then
	  write(0,*) "Usage : intratorsion2 <HISTORYfile> <OUTPUTfile> <sp> <i1> <j1> <k1> <l1> <i2> <j2> <k2> <l2>"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(i4)") sp
	call getarg(4,temp); read(temp,"(i4)") t1(1)
	call getarg(5,temp); read(temp,"(i4)") t1(2)
	call getarg(6,temp); read(temp,"(i4)") t1(3)
	call getarg(7,temp); read(temp,"(i4)") t1(4)
	call getarg(8,temp); read(temp,"(i4)") t2(1)
	call getarg(9,temp); read(temp,"(i4)") t2(2)
	call getarg(10,temp); read(temp,"(i4)") t2(3)
	call getarg(11,temp); read(temp,"(i4)") t2(4)

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	write(0,"(A,I5)") "Target species is ",sp
	write(0,"(a,4(i3),a,4(i3))") "Calculating torsion angle map between ",t1(1),t1(2),t1(3),t1(4), " and ", t2(1),t2(2),t2(3),t2(4)
	
	allocate(ijkl(-ntorsionbins:ntorsionbins,-ntorsionbins:ntorsionbins),stat=status); if (status.GT.0) stop "Allocation error for ijkl()"

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
	
	aoff1 = s_start(sp)-1
	do s1 = 1,s_nmols(sp)     ! Loop over all molecules in species

	  psi1 = torsion(aoff1+t1(1),aoff1+t1(2),aoff1+t1(3),aoff1+t1(4))
	  psi2 = torsion(aoff1+t2(1),aoff1+t2(2),aoff1+t2(3),aoff1+t2(4))

	  if (psi1.lt.0.0) then
	    bin1 = int(psi1/torsionbin) - 1
	    if (bin1.lt.-ntorsionbins) bin1 = -ntorsionbins
	  else
	    bin1 = int(psi1/torsionbin)
	    if (bin1.ge.ntorsionbins) bin1 = ntorsionbins-1
	  end if
	  if (psi2.lt.0.0) then
	    bin2 = int(psi2/torsionbin) - 1
	    if (bin2.lt.-ntorsionbins) bin2 = -ntorsionbins
	  else
	    bin2 = int(psi2/torsionbin)
	    if (bin2.ge.ntorsionbins) bin2 = ntorsionbins-1
	  end if
	  ijkl(bin1,bin2) = ijkl(bin1,bin2) + 1

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
	ijkl = ijkl / (nframes * s_nmols(sp) * s_nmols(sp))

	a(1) = char(48+t1(1)/10)//char(48+MOD(t1(1),10))
	a(2) = char(48+t1(2)/10)//char(48+MOD(t1(2),10))
	a(3) = char(48+t1(3)/10)//char(48+MOD(t1(3),10))
	a(4) = char(48+t1(4)/10)//char(48+MOD(t1(4),10))
	a(5) = char(48+t2(1)/10)//char(48+MOD(t2(1),10))
	a(6) = char(48+t2(2)/10)//char(48+MOD(t2(2),10))
	a(7) = char(48+t2(3)/10)//char(48+MOD(t2(3),10))
	a(8) = char(48+t2(4)/10)//char(48+MOD(t2(4),10))

	! Write torsion histogram ijkl
	resfile=basename(1:baselen)//a(1)//"-"//a(2)//"-"//a(3)//"-"//a(4)//"_"//a(5)//"-"//a(6)//"-"//a(7)//"-"//a(8)//".tors2"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=-ntorsionbins,ntorsionbins-1
	  write(9,"(f9.3,2x,f12.8)") ((n+0.5)*torsionbin,ijkl(n,m),m=-ntorsionbins,ntorsionbins-1)
	  write(9,*) ""
	end do
	close(9)

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end program intratorsion2

	real*8 function torsion(a, b, c, d)
	use dlprw; use utility
	implicit none
	real*8, parameter :: pi = 3.14159265358979d0, radcon = 57.29577951d0
	integer :: a, b, c, d
	real*8 :: i(3), j(3), k(3), l(3), v(3), vecji(3), vecjk(3), veckj(3), veckl(3), xp1(3), xp2(3)
	real*8 :: norm, dp, rij, tx, ty, tz, ktx, kty, ktz, mag1, mag2
	real*8 :: dist, magnitude, dotproduct

	! Grab coordinates of the atoms
	i(1) = xpos(a)
	i(2) = ypos(a)
	i(3) = zpos(a)

	j(1) = xpos(b)
	j(2) = ypos(b)
	j(3) = zpos(b)

	k(1) = xpos(c)
	k(2) = ypos(c)
	k(3) = zpos(c)

	l(1) = xpos(d)
	l(2) = ypos(d)
	l(3) = zpos(d)

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
	if (dp.lt.0) torsion = -torsion
	end function torsion

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

