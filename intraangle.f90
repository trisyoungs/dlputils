!	** intraangle **
!	Calculate angle between vectors on the same molecule

	program intraangle
	use dlprw; use utility
	implicit none
	real*8, parameter :: pi = 3.14159265358979d0
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*2 :: a1, a2, a3
	character*20 :: temp
	integer :: status,  nargs, success, baselen
	integer :: nanglebins,aoff1,n,m,s1,bin,bin2,nframes,numadded,sp,i,j,k
	integer :: iargc
	real*8 :: ri(3), rj(3), rk(3), rji(3), rjk(3), v(3),dp, anglebin, norm, angle, binsum
	real*8, allocatable :: ijk(:)

	anglebin=0.1	! In degrees

	nargs = iargc()
	if (nargs.ne.6) then
	  write(0,*) "Usage : intraangle <HISTORYfile> <OUTPUTfile> <sp> <i> <j> <k>"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(i4)") sp
	call getarg(4,temp); read(temp,"(i4)") i
	call getarg(5,temp); read(temp,"(i4)") j
	call getarg(6,temp); read(temp,"(i4)") k

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	nanglebins = 180.0 / anglebin
	write(0,"(A,F6.3,A)") "Using binwidths of ",anglebin,"Degrees"
	write(0,"(A,I5,A)") "There will be ",nanglebins," angle bins."
	write(0,"(A,I5)") "Target species is ",sp
	write(0,"(a,3i3)") "Calculating angle between vectors from atoms ",i,j,k
	
	allocate(ijk(nanglebins),stat=status); if (status.GT.0) stop "Allocation error for ijk()"

	! Initialise the arrays...
	ijk = 0.0

	! XXXX
	! XXXX Main routine....
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

	  ! Grab coordinates of the first vector (ri->rj)
	  ri(1) = xpos(aoff1+i-1)
	  ri(2) = ypos(aoff1+i-1)
	  ri(3) = zpos(aoff1+i-1)

	  rj(1) = xpos(aoff1+j-1)
	  rj(2) = ypos(aoff1+j-1)
	  rj(3) = zpos(aoff1+j-1)

	  rk(1) = xpos(aoff1+k-1)
	  rk(2) = ypos(aoff1+k-1)
	  rk(3) = zpos(aoff1+k-1)

	  ! PBC atoms to j
	  call pbc(ri(1),ri(2),ri(3),rj(1),rj(2),rj(3),ri(1),ri(2),ri(3))
	  call pbc(rk(1),rk(2),rk(3),rj(1),rj(2),rj(3),rk(1),rk(2),rk(3))

	  ! Calculate vectors and normalise
	  rji = ri - rj
	  rji = rji / sqrt(sum(rji*rji))
	  rjk = rk - rj
	  rjk = rjk / sqrt(sum(rjk*rjk))

	  ! dist = sqrt(irl(1)*irl(1) + irl(2)*irl(2) + irl(3)*irl(3))

	  !
	  ! Analysis
	  !
	  ! Angle ri-rj-rl
	  dp = dot_product(rji, rjk)
	  angle = safeAngle(dp)
	  bin = int(angle * (1.0 / anglebin)) + 1
	  ijk(bin) = ijk(bin)+1


	  ! Global counter
	  numadded = numadded+1

	  aoff1 = aoff1 + s_natoms(sp)
	end do   ! End main loop over all atoms of species1.

	if (nframes.EQ.1) write(0,*) "numadded:=",numadded
	if (nframes.EQ.1) then
	  write(0,"(A,I2,A,I4,A)") "PRDF of atoms about ",sp," : averaged over ",s_nmols(sp)," molecules."
	end if

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

	! Normalise to number of frames
	ijk = ijk / nframes

	binsum = sum(ijk)

	a1 = char(48+i/10)//char(48+MOD(i,10))
	a2 = char(48+j/10)//char(48+MOD(j,10))
	a3 = char(48+k/10)//char(48+MOD(k,10))

	! Write histogram ijk
	resfile=basename(1:baselen)//a1//"-"//a2//"-"//a3//".ijk"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,nanglebins
	  write(9,"(F10.3,2(3x,e12.6))") anglebin*(n-0.5),ijk(n),ijk(n)/binsum
	end do
	close(9)

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end program intraangle

