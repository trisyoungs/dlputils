!	** cdf2 **
!	Calculate cylindrical distribution functions between the centres-of-mass of
!	molecules and a specified vector.
!	Calculates the density at each frame of the simulation (for NPT simulations)
!	Performs moving slice average function.
!	NOTE: Only valid / meaningful when pore is along z axis of cell.

	program cyldf2
	use dlprw; use utility; use parse
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	logical :: altheader = .FALSE.
	integer :: nbins,n,m,a1,s1,m1,baselen,bin,nframes,success,nargs,numadded,framestodo = -1,framestodiscard = 0,framesdone, compairs(10,2)
	integer :: nsamples = 100
	integer :: iargc
	real*8 :: dist,pos(3),integral,origin(3),vector(3),t1(3),t2(3),z,zmin,zmax
	real*8 :: binwidth, xyyx, xzzx, yzzy, denom, numdens, shellvol, slicewidth = -1.0, sampledelta
	real*8, allocatable :: cdf(:,:)

	binwidth=0.1   ! In Angstroms
	compairs = 0
	nargs = iargc()
	if (nargs.LT.4) stop "Usage : cdf2 <HISTORYfile> <OUTPUTfile> <ox> <oy> [-bin width] [-slice width] [-nsamples n] [-header hisfile] [-frames n] [-discard n] [-compair sp i j]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(f20.14)") origin(1)
	call getarg(4,temp); read(temp,"(f20.14)") origin(2)
	origin(3) = 0.0
	vector(1) = 0.0
	vector(2) = 0.0
	vector(3) = 1.0
	write(0,"(a,3f10.4)") "Vector origin is : ", origin

	! Check line parameters
	denom = dot_product(vector,vector)
	if (denom.lt.1.0e-6) stop "Invalid line vector supplied."
	
	n = 4
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-bin")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") binwidth
	      write(0,"(A,f6.2)") "Binwidth set to ",binwidth
	    case ("-header")
	      n = n + 1; call getarg(n,altheaderfile)
	      write(0,"(A,I4)") "Alternative header file supplied."
	      altheader = .TRUE.
	    case ("-frames")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
	      write(0,"(A,I4)") "Frames to process: ",framestodo
	    case ("-discard")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodiscard
	      write(0,"(A,I4)") "Frames to discard at start: ",framestodiscard
	    case ("-compair")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") s1
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(s1,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(s1,2)
	      write(0,"(A,3I4)") "Using COMpair for species ",s1, compairs(s1,:)
	    case ("-slice")
	      n = n + 1; slicewidth = getargr(n)
	      write(0,"(A,f6.2)") "Slice width set to ",slicewidth
	    case ("-nsamples")
	      n = n + 1; nsamples = getargi(n)
	      write(0,"(A,f6.2)") "Number of samples set to ",nsamples
	    case default
	      write(0,"(a,a)") "Unrecognised command line option:",temp
	      stop
	  end select
	end do

	! Open and check the files...
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798

	! Now, read in the history header so that we have cell()
	! If this fails then we may have a restarted trajectory. Continue, but only if
	! a header can be read in from the specified alternative history file..
	call openhis(hisfile,10)
	if (readheader().EQ.-1) then
	  if (altheader) then
	    write(0,*) "Restarted trajectory:"
	    close(dlpun_his)
	    call openhis(altheaderfile,10)
	    if (readheader().EQ.-1) goto 797
	    close(dlpun_his)
	    call openhis(hisfile,10)
	  else
	    goto 797
	  end if
	end if

	! If slice width is less than zero, set it to one tenth of the vector
	if (slicewidth < 0.0) slicewidth = 0.1*(vector(1)*cell(1) + vector(2)*cell(5) + vector(3)*cell(9))
	sampledelta = ((vector(1)*cell(1) + vector(2)*cell(5) + vector(3)*cell(9)) - slicewidth) / nsamples
	write(0,"(A,I5,A,F6.3,A)") "There will be ",nsamples," samples of slices of width ",slicewidth," Angstroms."
	write(0,"(A,F6.3)") "Sampling delta is ",sampledelta

	nbins = maxval(cell) / binwidth + 1
	write(0,"(A,I5,A,F6.3,A)") "There will be ",nbins," histogram bins of ",binwidth," Angstroms."
	allocate(cdf(nspecies,nbins))
	cdf = 0.0

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
100	nframes=0
	framesdone = 0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes, framesdone
	if (nframes.le.framestodiscard) goto 101

	framesdone = framesdone + 1

	call calc_com

	do s1=1,nspecies
	  ! If compairs were specified, use that instead of COM
	  if (compairs(s1,1).ne.0) then
	    do m1=1,s_nmols(s1)
	      t1(1) = xpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,1)-1)
	      t1(2) = ypos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,1)-1)
	      t1(3) = zpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,1)-1)
	      t2(1) = xpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,2)-1)
	      t2(2) = ypos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,2)-1)
	      t2(3) = zpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,2)-1)
	      call pbc(t2(1),t2(2),t2(3),t1(1),t1(2),t1(3),pos(1),pos(2),pos(3))
	      comx(s1,m1) = (t1(1)+pos(1))*0.5
	      comy(s1,m1) = (t1(2)+pos(2))*0.5
	      comz(s1,m1) = (t1(3)+pos(3))*0.5
	    end do
	  end if

	  zmin = 0.0
	  zmax = slicewidth

	  do n=1,nsamples
	    do m1=1,s_nmols(s1)     ! Loop over all molecules of species...
	      ! Check z position, folding into range 0 - cell(9) is necessary
	      z = comz(s1,m1)
	      if (z.lt.0.0) z = z + cell(9)
	      if (z.gt.cell(9)) z = z - cell(9)
	      if (z.lt.zmin) cycle
	      if (z.gt.zmax) cycle

	      pos(1) = comx(s1,m1) - origin(1)
	      pos(2) = comy(s1,m1) - origin(2)
	      pos(3) = comz(s1,m1) - origin(3)
	      xyyx = vector(1)*pos(2) - vector(2)*pos(1)
	      xzzx = vector(1)*pos(3) - vector(3)*pos(1)
	      yzzy = vector(2)*pos(3) - vector(3)*pos(2)
	      t1(1) = vector(2)*xyyx + vector(3)*xzzx
	      t1(2) = vector(3)*yzzy - vector(1)*xyyx
	      t1(3) = -vector(1)*xzzx - vector(2)*yzzy
	      dist = sqrt(sum(t1*t1)) / denom
	      ! 'Add' this distance...
	      bin=INT(dist/binwidth)+1
	      if (bin.lt.nbins) cdf(s1,bin) = cdf(s1,bin) + 1.0
	    end do

	  ! Shift sample position
	  zmin = zmin + sampledelta
	  zmax = zmax + sampledelta

	  end do
	end do

	if (framesdone.eq.framestodo) goto 801

	! Next frame
	goto 101

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

797	write(0,*) "No header found in history file. If a restarted trajectory, use '-header'"
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
	  if (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="rdfresults."
	  baselen=11
	else
	  basename=hisfile(1:baselen)
	endif

	! Normalise w.r.t. number of frames and number of samples
	cdf = cdf / framesdone
	cdf = cdf / nsamples

	do s1=1,nspecies
	  resfile=basename(1:baselen)//"cdf2"//CHAR(48+s1)
	  open(unit=9,file=resfile,form="formatted")
	  write(9,"('# Origin / Vector = ',6f10.4)") origin, vector
	  if (compairs(s1,1).eq.0) then
	    write(9,"(a,i2,a)") "# Species ",s1," calculated using all atoms for COM."
	  else
	    write(9,"(a,i2,a,i3,i3)") "# Species ",s1," calculated using two atoms for COM: ", compairs(s1,1), compairs(s1,2)
	  end if
	  integral = 0.0
	  ! Normalise the RDFs with respect to the number of frames, cylindrical shell volume, and number density of species
	  numdens = s_nmols(s1) / volume(cell)
	  do n=1,nbins
	    integral = integral + cdf(s1,n)
	    ! Note - cylinder shell volume will only be correct for line vector coincident with cell vectors
	    !shellvol = PI*(n*binwidth)**2 - PI*((n-1)*binwidth)**2
	    shellvol = 3.141592654d0*binwidth*binwidth*( n*n - (n-1)*(n-1) ) * slicewidth
	    cdf(s1,n) = cdf(s1,n) / (shellvol * numdens)
	    ! write(9,"(F6.3,3x,F12.8)") (n*binwidth)-binwidth/2.0,nrdf(a1,n)
	    write(9,"(F10.4,2(3x,F12.8))") (n-0.5)*binwidth, cdf(s1,n), integral
	  end do
	  close(9)
	end do
	write(0,*) "Finished."
999	close(10)
	close(13)
	end program cyldf2

