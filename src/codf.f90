!	** cdf **
!	Calculate cylindrical, orientational distribution functions between a specified vector
!	and the axes of the specified molecule

	program codf
	use dlprw; use utility; use parse; use IList
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	real*8, parameter :: pi = 3.14159265358979d0
	real*8, parameter :: radcon = 57.29577951d0
	logical :: altheader = .FALSE.
	integer :: i,n,m,a1,sp,m1,baselen,bin,nframes,success,nargs,numadded,framestodo = -1,framestodiscard = 0,framesdone, compairs(10,2)
	integer :: nDistBins, nAngleBins, distBin, angleBin, axis
	integer :: iargc
	real*8 :: dist,pos(3),origin(3),vector(3),t1(3),t2(3),dp,angle,norm
	real*8 :: distBinWidth, angleBinWidth, xyyx, xzzx, yzzy, denom, numdens, shellvol
	type(IntegerList) :: targetsp
	real*8, allocatable :: cdf(:,:,:), cdfangles(:,:,:,:)

	distBinWidth=0.1   ! In Angstroms
	angleBinWidth = 1.0
	compairs = 0
	nargs = iargc()
	if (nargs.LT.8) stop "Usage : codf <HISTORYfile> <OUTPUTfile> <ox> <oy> <oz> <vx> <vy> <vz> <targetsp> [-bin width] [-anglebin width] [-header hisfile] [-frames n] [-discard n] [-axis sp x1 x2 y1 y2] [-compair i j]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	origin(1) = getargr(3)
	origin(2) = getargr(4)
	origin(3) = getargr(5)
	vector(1) = getargr(6)
	vector(2) = getargr(7)
	vector(3) = getargr(8)
	call getarg(9,temp); if (.not.parseIntegerList(temp, targetsp)) stop "Failed to parse targetsp list."
	write(0,"(a,i4)") "Target species are : ", targetsp%items(1:targetsp%n)
	write(0,"(a,3f10.4)") "Vector origin is : ", origin
	write(0,"(a,3f10.4)") "Vector is : ", vector

	! Check line parameters
	denom = dot_product(vector,vector)
	if (denom.lt.1.0e-6) stop "Invalid line vector supplied."

	! Open and check the output file so we know the number of species etc.
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798

	call alloc_axis()
	
	n = 9
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-anglebin")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") angleBinWidth
	      write(0,"(A,f6.2)") "Angle binwidth set to ", angleBinWidth
	    case ("-axis") 
	      n = n + 1; sp = getargi(n)
	      n = n + 1; axesAatoms(sp,1) = getargi(n)
	      n = n + 1; axesAatoms(sp,2) = getargi(n)
	      n = n + 1; axesAatoms(sp,3) = getargi(n)
	      n = n + 1; axesAatoms(sp,4) = getargi(n)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp," calculated from: X=",axesAatoms(sp,1),"->", &
	        & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
	      do i=1,3
		if ((axesAatoms(sp,i).lt.1).or.(axesAatoms(sp,i).gt.s_natoms(sp))) stop "Atom id out of range for axes on this species!"
	      end do
	      axesAdefined(sp) = .true.
	    case ("-bin")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") distBinWidth
	      write(0,"(A,f6.2)") "Distance binwidth set to ", distBinWidth
	    case ("-compair")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") sp
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(sp,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(sp,2)
	      write(0,"(A,3I4)") "Using COMpair for species ",sp, compairs(sp,:)
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
	    case default
	      write(0,"(a,a)") "Unrecognised command line option:",temp
	      stop
	  end select
	end do

	! Do we have all relevant axes definitions?
	do n=1,targetsp%n
	  sp = targetsp%items(n)
	  ! Check that the necessary molecules have had their axes defined
	  if (axesAdefined(sp)) then
	    write(6,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp," calculated from: X=",axesAatoms(sp,1),"->", &
	      & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
	  else
	    stop "Axes must be defined on all target species."
	  end if
	end do

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

	nAngleBins = 180.0 / angleBinWidth + 1
	nDistBins = maxval(cell) / distBinWidth + 1
	write(0,"(A,I5,A,F6.3,A)") "There will be ",nDistBins," distance histogram bins of ", distBinWidth," Angstroms."
	write(0,"(A,I5,A,F6.3,A)") "There will be ",nAngleBins," angle histogram bins of ", angleBinWidth," degrees."
	allocate(cdfangles(targetsp%n,3,nDistBins,nAngleBins))
	allocate(cdf(targetsp%n,nDistBins,0:6))
	cdf = 0.0
	cdfangles = 0.0

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

	! Generate centre of mass and all molecular axes
	call calc_com()
	call genaxes()

	do n=1,targetsp%n
	  sp = targetsp%items(n)
	  ! If compairs were specified, use that instead of COM
	  if (compairs(sp,1).ne.0) then
	    do m1=1,s_nmols(sp)
	      t1(1) = xpos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,1)-1)
	      t1(2) = ypos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,1)-1)
	      t1(3) = zpos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,1)-1)
	      t2(1) = xpos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,2)-1)
	      t2(2) = ypos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,2)-1)
	      t2(3) = zpos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,2)-1)
	      call pbc(t2(1),t2(2),t2(3),t1(1),t1(2),t1(3),pos(1),pos(2),pos(3))
	      comx(sp,m1) = (t1(1)+pos(1))*0.5
	      comy(sp,m1) = (t1(2)+pos(2))*0.5
	      comz(sp,m1) = (t1(3)+pos(3))*0.5
	    end do
	  end if

	  do m1=1,s_nmols(sp)     ! Loop over all molecules of the target species

	    ! Calculate the distance from the reference vector
	    pos(1) = comx(sp,m1) - origin(1)
	    pos(2) = comy(sp,m1) - origin(2)
	    pos(3) = comz(sp,m1) - origin(3)
	    xyyx = vector(1)*pos(2) - vector(2)*pos(1)
	    xzzx = vector(1)*pos(3) - vector(3)*pos(1)
	    yzzy = vector(2)*pos(3) - vector(3)*pos(2)
	    t1(1) = vector(2)*xyyx + vector(3)*xzzx
	    t1(2) = vector(3)*yzzy - vector(1)*xyyx
	    t1(3) = -vector(1)*xzzx - vector(2)*yzzy
	    dist = sqrt(sum(t1*t1)) / denom
	    t1 = t1 / sqrt(sum(t1*t1))

	    ! Accumulate vector dot-products with distance vector (perpendicular to reference vector)
	    ! Calculate angles between distance vector (perpendicular to reference vector) and molecule axes
	    distBin=INT(dist/distBinWidth)+1
	    if (distBin.lt.nDistBins) then
	      cdf(n,distBin,0) = cdf(n,distBin,0) + 1.0

	      ! -- X
	      dp = vec3DotProduct(t1, axesA(sp,m1,1:3))
	      angle = safeAngle(dp)
	      angleBin = INT( angle / angleBinWidth) + 1
	!write(0,*) sp, m1, dist, dp, angleBin
	      cdf(n,distBin,1) = cdf(n,distBin,1) + dp
	      cdf(n,distBin,4) = cdf(n,distBin,4) + dabs(dp)
	      cdfangles(n,1,distBin,angleBin) = cdfangles(n,1,distBin,angleBin) + 1.0

	      ! -- Y
	      dp = vec3DotProduct(t1, axesA(sp,m1,4:6))
	      angle = safeAngle(dp)
	      angleBin = INT( angle / angleBinWidth) + 1
	      cdf(n,distBin,2) = cdf(n,distBin,2) + dp
	      cdf(n,distBin,5) = cdf(n,distBin,5) + dabs(dp)
	      cdfangles(n,2,distBin,angleBin) = cdfangles(n,2,distBin,angleBin) + 1.0

	      ! -- Z
	      dp = vec3DotProduct(t1, axesA(sp,m1,7:9))
	      angle = safeAngle(dp)
	      angleBin = INT( angle / angleBinWidth) + 1
	      cdf(n,distBin,3) = cdf(n,distBin,3) + dp
	      cdf(n,distBin,6) = cdf(n,distBin,6) + dabs(dp)
	      cdfangles(n,3,distBin,angleBin) = cdfangles(n,3,distBin,angleBin) + 1.0
	    end if
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

	do n=1,targetsp%n
	  sp = targetsp%items(n)

	  ! Dot product histograms
	  resfile=basename(1:baselen)//"codf"//CHAR(48+sp)
	  open(unit=9,file=resfile,form="formatted")
	  write(9,"('# Origin / Vector = ',6f10.4)") origin, vector
	  write(9,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "# Local axes for species ",sp," calculated from: X=",axesAatoms(sp,1),"->", &
	    & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
	  if (compairs(sp,1).eq.0) then
	    write(9,"(a,i2,a)") "# Species ",sp," calculated using all atoms for COM."
	  else
	    write(9,"(a,i2,a,i3,i3)") "# Species ",sp," calculated using two atoms for COM: ", compairs(sp,1), compairs(sp,2)
	  end if

	  do i=1,nDistBins
	    if (cdf(n,i,0).gt.0.5) cdf(n,i,1:6) = cdf(n,i,1:6) / cdf(n,i,0)
	    write(9,"(F10.4,6(3x,F12.8),3x,e12.5)") (i-0.5)*distBinWidth, (cdf(n,i,m), m=1,6), cdf(n,i,0)
	  end do
	  close(9)

	  do axis=1,3

	    ! Normalised histograms (sum to 1.0)

	    resfile=basename(1:baselen)//"codf"//CHAR(48+sp)//CHAR(119+axis)
	    open(unit=9,file=resfile,form="formatted")
	    write(9,"('# Origin / Vector = ',6f10.4)") origin, vector
	    write(9,"('# RDelta / RNBins = ',f10.4,i5)") distBinWidth, nDistBins
	    write(9,"('# AngDelta / AngNBins = ',f10.4,i5)") angleBinWidth, nAngleBins
	    write(9,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "# Local axes for species ",sp," calculated from: X=",axesAatoms(sp,1),"->", &
	      & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
	    if (compairs(sp,1).eq.0) then
	      write(9,"(a,i2,a)") "# Species ",sp," calculated using all atoms for COM."
	    else
	      write(9,"(a,i2,a,i3,i3)") "# Species ",sp," calculated using two atoms for COM: ", compairs(sp,1), compairs(sp,2)
	    end if

	    do i=1,nDistBins
	      ! Normalise histograms (per distance)
	      norm = sum(cdfangles(n,axis,i,:))
	      if (norm > 0.5) cdfangles(n,axis,i,:) = cdfangles(n,axis,i,:) / norm
	    
	      do m=1,nAngleBins
	        write(9,"(F10.4,3x,F12.8)") (m-0.5)*angleBinWidth, cdfangles(n,axis,i,m) / dsin(((m-0.5)*angleBinWidth)/RADCON)
	      end do
	      write(9,*) ""
	    end do
	    close(9)

	    ! Histograms normalised to cylindrical shell volume
	    ! Assumes length of cylinder is same as unit cell C length

	    resfile=basename(1:baselen)//"codf_cylnorm"//CHAR(48+sp)//CHAR(119+axis)
	    open(unit=9,file=resfile,form="formatted")
	    write(9,"('# Origin / Vector = ',6f10.4)") origin, vector
	    write(9,"('# RDelta / RNBins = ',f10.4,i5)") distBinWidth, nDistBins
	    write(9,"('# AngDelta / AngNBins = ',f10.4,i5)") angleBinWidth, nAngleBins
	    write(9,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "# Local axes for species ",sp," calculated from: X=",axesAatoms(sp,1),"->", &
	      & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
	    if (compairs(sp,1).eq.0) then
	      write(9,"(a,i2,a)") "# Species ",sp," calculated using all atoms for COM."
	    else
	      write(9,"(a,i2,a,i3,i3)") "# Species ",sp," calculated using two atoms for COM: ", compairs(sp,1), compairs(sp,2)
	    end if

	    do i=1,nDistBins
	      ! Normalise to cylindrical shell volume, assuming unit cell C is cylinder length
	      norm = (3.141592654d0*distBinWidth*distBinWidth*( n*n - (n-1)*(n-1) ) * cell(9)) * (s_nmols(sp) / volume(cell))
	      cdfangles(n,axis,i,:) = (cdfangles(n,axis,i,:) / framesdone) / norm
	    
	      do m=1,nAngleBins
	        write(9,"(F10.4,3x,F12.8)") (m-0.5)*angleBinWidth, cdfangles(n,axis,i,m) / dsin(((m-0.5)*angleBinWidth)/RADCON)
	      end do
	      write(9,*) ""
	    end do
	    close(9)

	  end do
	
	end do

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program codf

