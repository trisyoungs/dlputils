!	** intradist **
!	Calculate angle between vectors on the same molecule

	program intradist
	use dlprw; use utility
	implicit none
	real*8, parameter :: pi = 3.14159265358979d0
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*2 :: a1, a2, a3
	character*20 :: temp
	integer :: status,  nargs, success, baselen
	integer :: nbins,aoff1,n,s1,bin,nframes,numadded,sp,i,j
	integer :: iargc
	real*8 :: ri(3), rj(3), rij(3), binwidth, binsum
	real*8, allocatable :: ij(:)

	binwidth=0.01

	nargs = iargc()
	if (nargs.ne.5) then
	  write(0,*) "Usage : intradist <HISTORYfile> <OUTPUTfile> <sp> <i> <j>"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(i4)") sp
	call getarg(4,temp); read(temp,"(i4)") i
	call getarg(5,temp); read(temp,"(i4)") j

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	nbins = 3.0 / binwidth
	write(0,"(A,F6.3)") "Using binwidths of ",binwidth
	write(0,"(A,I5,A)") "There will be ",nbins," distance bins."
	write(0,"(A,I5)") "Target species is ",sp
	write(0,"(a,2i3)") "Calculating distance between atoms ",i,j
	
	allocate(ij(nbins),stat=status); if (status.GT.0) stop "Allocation error for ij()"

	! Initialise the arrays...
	ij = 0.0

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

	  ! Grab coordinates of the atoms
	  ri(1) = xpos(aoff1+i-1)
	  ri(2) = ypos(aoff1+i-1)
	  ri(3) = zpos(aoff1+i-1)

	  rj(1) = xpos(aoff1+j-1)
	  rj(2) = ypos(aoff1+j-1)
	  rj(3) = zpos(aoff1+j-1)

	  ! PBC atom i to j
	  call pbc(ri(1),ri(2),ri(3),rj(1),rj(2),rj(3),ri(1),ri(2),ri(3))

	  ! Calculate vectors and normalise
	  rij = ri - rj

	  !
	  ! Analysis
	  !
	  bin = int(sqrt(sum(rij*rij)) * (1.0 / binwidth)) + 1
	  ij(bin) = ij(bin)+1

	  ! Global counter
	  numadded = numadded+1

	  aoff1 = aoff1 + s_natoms(sp)
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

	! Normalise to number of frames
	ij = ij / nframes

	a1 = char(48+i/10)//char(48+MOD(i,10))
	a2 = char(48+j/10)//char(48+MOD(j,10))

	binsum = sum(ij)

	! Write histogram ij
	resfile=basename(1:baselen)//a1//"-"//a2//".ij"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,nbins
	  write(9,"(F10.3,2(3x,e12.6))") binwidth*(n-0.5),ij(n),ij(n)/binsum
	end do
	close(9)

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end program intradist

