!	** staticsk **
!	Compute the static total structure factor S(k)

	program staticsk
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	logical :: npt = .FALSE.
	integer :: baselen,nframes,success,nargs,framestodo
	integer :: i,j,k,n, nbins, bin, found, kx, ky, kz, nvec
	integer, allocatable :: numadded(:)
	integer :: iargc
	real*8 :: binwidth,factor,kcut,magx,magy,magz,mag,kmag(20000)
	real*8, allocatable :: sk(:),rxxx(:),ryyy(:),rzzz(:)
	complex :: kvectors(20000)
	complex :: vecx, vecy, vecz
	real*8, parameter :: pi = 3.14159265358979d0

	binwidth=0.1
	kcut = 5.0    ! Reciprocal space cutoff (box integers)

	nargs = iargc()
	if (nargs.LT.2) stop "Usage : staticsk <HISTORYfile> <OUTPUTfile> [-bin binwidth] [-frames nframes] [-npt] [-kcut cutoff]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	if (nargs.GE.3) then
	  n = 3
	  do
	    call getarg(n,temp)
	    select case (temp)
	      case ("-bin"); n = n + 1; call getarg(n,temp); read(temp,"(F20.10)") binwidth
	      case ("-kcut"); n = n + 1; call getarg(n,temp); read(temp,"(F20.10)") kcut
	      case ("-frames"); n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
	      case ("-npt"); npt = .TRUE.
	    end select
	    n = n + 1
	    if (n.GT.nargs) exit
	  end do
	end if
	
	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell() and also the atom names
        if (readheader().EQ.-1) goto 799

	nbins = kcut / binwidth + 1
	allocate(sk(nbins))
	allocate(numadded(nbins))
	allocate(rxxx(natms))
	allocate(ryyy(natms))
	allocate(rzzz(natms))

	sk = 0.0
	numadded = 0.0

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
	
	! Calculate reciprocal cell for this configuration (if necessary)
	if (npt.OR.(nframes.EQ.1)) call calc_rcell

	! Convert atomic coordinates to reciprocal space coordinates
	rxxx = rcell(1)*xpos + rcell(2)*ypos + rcell(3)*zpos
	ryyy = rcell(4)*xpos + rcell(5)*ypos + rcell(6)*zpos
	rzzz = rcell(7)*xpos + rcell(8)*ypos + rcell(9)*zpos

	! Generate reciprocal space vectors
	kx = kcut * 2.0
	ky = kx
	kz = kx
	nvec = 0
	do i=-kx,kx
	  do j=-ky,ky
	    do k=-kz,kz
	      magx = i*rcell(1) + j*rcell(2) + k*rcell(3)
	      magy = i*rcell(4) + j*rcell(5) + k*rcell(6)
	      magz = i*rcell(7) + j*rcell(8) + k*rcell(9)
	      mag = sqrt(magx**2 + magy**2 + magz**2)
	      if ((mag.LT.kcut).AND.(mag.GT.0.01)) then
		nvec = nvec + 1
		kmag(nvec) = mag
		kvectors(nvec) = (0.0,0.0)
		do n=1,natms
		  vecx=cmplx(cos(-rxxx(n)),sin(-rxxx(n)))
		  vecy=cmplx(cos(-ryyy(n)),sin(-ryyy(n)))
		  vecz=cmplx(cos(-rzzz(n)),sin(-rzzz(n)))
		  ! Sum the atomic kvector contributions into kvectors()
		  kvectors(nvec) = kvectors(nvec) + (vecx**i) * (vecy**j) * (vecz**k)
		end do
	      end if
	    end do
	  end do
	end do

	if (nframes.EQ.1) write(0,"(I6,A)") nvec," k-vectors generated from first configuration..."

	! Sum kvectors into bins
	do n=1,nvec
	  bin = kmag(n) / binwidth + 1
	  sk(bin) = sk(bin) + kvectors(n)*conjg(kvectors(n))
	  numadded(bin) = numadded(bin) + 1
	end do

	write(0,*) "FRAME ",nframes, "SK(30) ",sk(30)

	! Next frame
	if (nframes.EQ.framestodo) goto 120
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
	  basename="sfact."
	  baselen=6
	ELSE
	  basename=hisfile(1:baselen)
	endif

	resfile=basename(1:baselen)//"sk"
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,nbins
	  if (numadded(n).EQ.0) then
	    sk(n) = 0.0
	  else
	    factor = nframes * numadded(n) * natms
	    sk(n) = sk(n) / factor
	  end if
	  write(9,"(2F15.7)") (n-0.5)*binwidth,sk(n)
	end do
	close(9)

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end program staticsk
