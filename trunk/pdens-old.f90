!	** 3dpd.f90 **
!	(based on old 3ddist_v2)
!	Program to calculate the 3-dimensional distribution of different species about
!	one species' centres.

	program pdens3d
	use dlprw; use utility
	implicit none
	real*8, allocatable :: pdens(:,:,:,:)
	real*8, allocatable :: avggeom(:,:,:)		! Average species coordinates
	real*8, allocatable :: pdensintra(:,:,:,:)	! Intramolecular species distributions
	integer, allocatable :: pdensn(:), pdensa(:)	! Total and 'binned' mols
	integer, allocatable :: molflags(:)		! Per-frame map of sp1 mols to include in averaging
	integer, allocatable :: spexp(:)		! Expected species numbers
	integer :: grid = 20, pgrid = 50	! Grid in each direction for 3ddist and pdens
	real*8 :: delta = 0.5, pdelta = 0.15	! Default grid spacings
	logical :: molmap = .FALSE.		! Whether we're using a file of mol flags
	integer :: maxframes,maxbins,patoms(2),nmapframes
	integer :: baselen, nargs, molcount,norm1,norm2
	character*80 :: hisfile,outfile,basename,resfile,temp,flagfile
	character*8 :: discard
	integer :: success,calctype,bin,X,numadded,comtype,n1,n2,n3
	integer :: iatm, n, discardn, nframes,k, l, m, o, p
	real*8 :: px,py,pz, summ
	real*8 :: boxvolume
	real*8 :: tx,ty,tz,discardr,c1x,c1y,c1z,c2x,c2y,c2z
	integer :: centresp,sp1,m1,sp2,m2,startf, endf,skip
	integer :: iargc

	write(0,*) "*** 3Dpd"

	nargs = iargc()
	if (nargs.LT.2) then
	  write(*,*) "Usage: 3dpd <HIStory file> <OUTput file> [-options]"
	  write(*,*) "        [-centre sp]            Set species sp to be the central one"
	  write(*,*) "        [-axis sp x1 x2 y1 y2]  Atoms to use for axis calculation in species sp"
	  write(*,*) "        [-grid npoints]         Grid points in each direction for prob. densities (default = 20)"
	  write(*,*) "        [-delta spacing]        Spacing between grid points (default = 0.5 Angstrom)"
	  write(*,*) "        [-pgrid npoints]        Grid points in each direction for species densities (default = 50)"
	  write(*,*) "        [-pdelta spacing]       Spacing between pgrid points (default = 0.15 Angstrom)"
	  write(*,*) "        [-start frameno]        Trajectory frame to start calculations (default = 1)"
	  write(*,*) "        [-end frameno]          Trajectory frame to end calculations on (default = last)"
	  write(*,*) "        [-molmap <file>]        Formatted file (I5,1X,n(I2,1x) of frame,mols to use of sp1"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
	end if
	n = 10
	write(0,"(A,A)") "History file : ",hisfile
	call openhis(hisfile,n)
	write(0,"(A,A)") " Output file : ",outfile
	if (outinfo(outfile,1).EQ.-1) goto 798

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  IF (hisfile(n:n).EQ.".") THEN
	    baselen=n
	    goto 50
	  endIF
	end do
50	if (baselen.EQ.-1) THEN
	  basename="3ddistresults."
	  baselen=14
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=15,file=basename(1:baselen)//"3dpdout",form="formatted",status="replace")

	! Set some variable defaults before we read in any command-line arguments
	allocate(aa(nspecies,4))
	aa = 0
	aa(1,1:4) = (/ 1,3,2,2 /)
	centresp = 1
	startf = 1
	endf = -1
	molcount = 0

	n = 2
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-axis") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp1
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp1,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp1,2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp1,3)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp1,4)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp1," calculated from: X=",aa(sp1,1),"->", &
	        & aa(sp1,2),", Y=0.5X->0.5(r(",aa(sp1,3),")->r(",aa(sp1,4),"))"
	    case ("-grid")
	      n = n + 1; call getarg(n,temp); read(temp,"(I4)") grid
	      write(0,"(A,I4)") "Grid points in each XYZ = ",grid
	    case ("-pgrid")
	      n = n + 1; call getarg(n,temp); read(temp,"(I4)") pgrid
	      write(0,"(A,I4)") "PGrid points in each XYZ = ",pgrid
	    case ("-delta")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") delta
	      write(0,"(A,F6.3)") "Grid spacing = ",delta
	    case ("-pdelta")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") pdelta
	      write(0,"(A,F6.3)") "PGrid spacing = ",pdelta
	    case ("-centre")
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") centresp
	      write(0,"(A,I4)") "Central species is = ",centresp
	    case ("-start")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") startf
	      write(0,"(A,I4)") "Starting frame = ",startf
	    case ("-end")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") endf
	      write(0,"(A,I4)") "End frame = ",endf
	    case ("-molmap")
	      n = n + 1; call getarg(n,flagfile)
	      write(0,"(A,A)") "Using flags for species 1 molecules from : ",flagfile
	      open(unit=20,file=flagfile,form="formatted",status="old")
	      molmap = .TRUE.
	    case default
	      write(0,*) "Unrecognised command line option : ",temp
	      stop
	  end select
	end do
	
	! Print out a summary of the control variables to the output file.
	write(15,"(A,A)") "Input file: ",hisfile
	write(15,"(A,I3)") "Molecular species in file : ",nspecies
	do sp1=1,nspecies
	  write(15,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Axis for species ",sp1," calculated from : X=",aa(sp1,1),"->", &
	    & aa(sp1,2),", Y=0.5X->0.5(",aa(sp1,3),"->",aa(sp1,4),")"
	end do
	write(15,"(A,I4)") "Grid (intermolecular) points in each XYZ = ",grid
	write(15,"(A,I4)") "PGrid (intramolecular) points in each XYZ = ",pgrid
	write(15,"(A,F6.3)") "Grid spacing = ",delta
	write(15,"(A,F6.3)") "PGrid spacing = ",pdelta
	write(15,"(A,I4)") "Central species = ",centresp
	write(15,"(A,I4,A,I4)") "Frame range = ",startf," to ",endf

	! Probability density arrays
	allocate(pdens(nspecies,-grid:grid,-grid:grid,-grid:grid))
	allocate(pdensa(nspecies)); allocate(pdensn(nspecies))
	allocate(pdensintra(nspecies,-pgrid:pgrid,-pgrid:pgrid,-pgrid:pgrid))
	allocate(spexp(nspecies))
	! Average species etc.
	allocate(avggeom(nspecies,maxatoms,3))
	if (molmap) allocate(molflags(s_nmols(centresp)))
	! Clear the arrays before we start...
	write(0,*) "Clearing initial arrays..."
	pdens(:,:,:,:) = reshape(  (/ (0.0,n=1,nspecies*(2*grid+1)**3) /) , (/ nspecies,2*grid+1,2*grid+1,2*grid+1 /) )
	pdensa(:) = (/ (0,n=1,nspecies) /)
	pdensn(:) = (/ (0,n=1,nspecies) /)
	pdensintra(:,:,:,:) = reshape(  (/ (0,n=1,nspecies*(2*pgrid+1)**3) /) , (/ nspecies,2*pgrid+1,2*pgrid+1,2*pgrid+1 /) )
	spexp(:) = (/ (0,n=1,nspecies) /)
	avggeom(:,:,:) = reshape( (/ (0,n=1,nspecies*maxatoms*3) /) , (/ nspecies,maxatoms,3 /) )

	! Read the header of the history file...
	IF (readheader().EQ.-1) goto 799

	! Clear some variables...
	pdensn(1:nspecies) = (/ (0,n=1,nspecies) /)
	pdensa(1:nspecies) = (/ (0,n=1,nspecies) /)
	pdens(1:nspecies,-grid:grid,-grid:grid,-grid:grid) = reshape( (/ (0.0,n=1,nspecies * &
	  & (grid*2+1)**3) /) , (/ nspecies,grid*2+1,grid*2+1,grid*2+1 /) )

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
	skip=0
100	nframes=0
	nmapframes=0
	! If we're using a mapfile, read it here
101	if (molmap) then
	  read(20,"(I5,1x,40(I2,1x))",end=118,err=117) m1,(molflags(m2),m2=1,s_nmols(centresp))
	  ! if (m1.NE.nframes) stop "Mapfile missed a frame number!"
	end if
102	success=readframe()
	IF (success.EQ.1) goto 120  ! End of file encountered....
	IF (success.EQ.-1) goto 119  ! File error....
	skip = skip + 1
	nframes=nframes+1
	if (MOD(nframes,100).EQ.0) then
	  if (molmap) then
	    write(0,"(I6,2x,'(',I6,',',I6,')')") nframes,nmapframes,molcount
	  else
	    write(0,"(i6)") nframes
	  end if
	end if
	if (skip.LT.startf) goto 102

	if ((molmap).and.(m1.ne.nframes)) goto 102
	if (molmap) nmapframes = nmapframes + 1

	! Calculate all geometric centres....
	call calc_com
	! Generate all molecular axes
	call genaxis

	sp1 = centresp
	do m1=1,s_nmols(sp1)      ! Loop over all molecules of species 1...
	  ! If we're using a mapfile, decide whether to include this molecule
	  if (molmap) then
	    if (molflags(m1).EQ.0) cycle 
	  end if
	  molcount = molcount + 1	! Normalisation counter
	  do sp2=1,nspecies
	    ! Now loop over all molecules of second species....
	    do m2=1,s_nmols(sp2)
	      ! * Calculations begin here *
	      ! ***************************
	      ! 1) Calculate the distribution of species sp2 about species sp1
	      call pbc(axisox(sp2,m2),axisoy(sp2,m2),axisoz(sp2,m2),axisox(sp1,m1), &
		& axisoy(sp1,m1), axisoz(sp1,m1),tx,ty,tz)
	      tx=tx-axisox(sp1,m1)
	      ty=ty-axisoy(sp1,m1)
	      tz=tz-axisoz(sp1,m1)
	      !write(40,*) comx(k,n),comy(k,n),comz(k,n),axisox(l,m),axisoy(l,m),axisoz(l,m)
	      !write(40,*) tx,ty,tz
	      ! Apply a transformation to rotate into the axis of the molecule l
	      px=tx*axisx(sp1,m1,1) + ty*axisx(sp1,m1,2) + tz*axisx(sp1,m1,3)
	      py=tx*axisy(sp1,m1,1) + ty*axisy(sp1,m1,2) + tz*axisy(sp1,m1,3)
	      pz=tx*axisz(sp1,m1,1) + ty*axisz(sp1,m1,2) + tz*axisz(sp1,m1,3)
	write(0,"(3(a3,3F10.4,4x))") "OLD",tx,ty,tz,"NEW",px,py,pz
	      ! Use n1,n2,n3 to store the integers for the pdens() array
	      n1=NINT(px/delta)
	      n2=NINT(py/delta)
	      n3=NINT(pz/delta)
	      ! If any of the n's are over grid, then only increase the counter
	      if (MAX(ABS(n1),MAX(ABS(n2),ABS(n3))).GT.grid) then
	        !write(0,*) "REFUSED DIST",l,m,k,n,SQRT(tx**2 + ty**2 + tz**2)
	        !write(0,*) "P",px,py,pz
	        !write(0,*) "PDIST",SQRT(px**2 + py**2 + pz**2)
	        !write(0,*) tx,ty,tz
	        !write(0,*) n1,n2,n3
	!write(0,"(3F10.4,3I4)") px,py,pz,n1,n2,n3
	        pdensn(sp2)=pdensn(sp2)+1      ! No position, just count it....
	      else
		if ((sp1.eq.sp2).and.(m1.eq.m2)) then
		  ! Same molecule, so don't add at all
		else
	          ! Check to make sure the same molecule isn't being consider with itself  XXXXX
	          pdens(sp2,n1,n2,n3)=pdens(sp2,n1,n2,n3)+1
	          pdensn(sp2)=pdensn(sp2)+1
	          pdensa(sp2)=pdensa(sp2)+1
		end if
	      end if
	      p=p+s_natoms(sp2)
	      ! * Calculations end *
	      ! ********************
	    end do
	  end do
	end do

	! Calculate average molecules and intramolecular distributions.
	do sp1=1,nspecies
	  p=s_start(sp1)
	  do m1=1,s_nmols(sp1)
	    ! If we're using a mapfile, decide whether to include this molecule
	    if ((sp1.EQ.centresp).AND.(molmap)) then
	      if (molflags(m1).EQ.0) cycle
	    end if
	    do n=1,s_natoms(sp1)
	      ! PBC the atom
	      call pbc(xpos(p+n-1),ypos(p+n-1),zpos(p+n-1),axisox(sp1,m1), &
	  	& axisoy(sp1,m1),axisoz(sp1,m1),tx,ty,tz)
	      ! 'Zero' its position with respect to the axis centre...
	      tx=tx-axisox(sp1,m1)
	      ty=ty-axisoy(sp1,m1)
	      tz=tz-axisoz(sp1,m1)
	      ! Transform the coordinate into the local coordinate system...
	      px=tx*axisx(sp1,m1,1) + ty*axisx(sp1,m1,2) + tz*axisx(sp1,m1,3)
	      py=tx*axisy(sp1,m1,1) + ty*axisy(sp1,m1,2) + tz*axisy(sp1,m1,3)
	      pz=tx*axisz(sp1,m1,1) + ty*axisz(sp1,m1,2) + tz*axisz(sp1,m1,3)
	      ! Accumulate the position....
	      avggeom(sp1,n,1)=avggeom(sp1,n,1)+px
	      avggeom(sp1,n,2)=avggeom(sp1,n,2)+py
	      avggeom(sp1,n,3)=avggeom(sp1,n,3)+pz
	      ! Intramolecular distribution.....
	      n1=NINT(px/pdelta)
	      n2=NINT(py/pdelta)
	      n3=NINT(pz/pdelta)
	      pdensintra(sp1,n1,n2,n3)=pdensintra(sp1,n1,n2,n3)+1
	    end do
	    p=p+s_natoms(sp1)
	  end do
	end do

	! Next frame
	if ((endf.NE.-1).AND.(nframes+startf.GT.endf)) goto 120
	goto 101

117	write(0,*) "Error reading mapfile."
	write(0,*) "Selected ",nmapframes," from trajectory before error."
	goto 120
118	write(0,*) "Reached end of mapfile."
	write(0,*) "Selected ",nmapframes," from trajectory."
	goto 120
119	write(0,*) "HISTORY file ended prematurely!"
	write(0,"(A,I5,A)") "Managed ",nframes," frames before error."
	write(15,"(A)") "HISTORY file ended prematurely!"
	write(15,"(A,I5,A)") "Managed ",nframes," frames before error."
120	write(0,*) "Finished."
	write(15,"(A)") "Finished."
	write(15,"(A,I5,A)") "Averages will be taken over ",nframes," frames."

	! ### Normalise the data

	! Set normalisation over frames and central species molecules
	do sp1=1,nspecies

	  if (molmap) then
	    norm1 = molcount
	    if (sp1.eq.centresp) then
	      norm2 = molcount
	    else
	      norm2 = nframes * s_nmols(sp1)
	    end if
	  else
	    norm1 = nframes * s_nmols(centresp)
	    norm2 = nframes * s_nmols(sp1)
	  end if
	  write(0,"(a,i1,a,i7,a,i7,a)") "Normalisation factors for species ",sp1," are: ",norm1," (external) and ",norm2," (internal)"

	  ! Calculate expected species numbers
	  spexp(sp1) = s_nmols(centresp) * s_nmols(sp1) * nframes
	  if (centresp.eq.sp1) spexp(sp1) = spexp(sp1) - s_nmols(sp1) * nframes

	  ! Species density about central species
	  do n1=-grid,grid
	    do n2=-grid,grid
	      do n3=-grid,grid
	        ! pdens(sp1,n1,n2,n3)=pdens(sp1,n1,n2,n3)/nframes/s_nmols(centresp)/(delta**3)
	        pdens(sp1,n1,n2,n3)=pdens(sp1,n1,n2,n3)/norm1/(delta**3)
	      end do
	    end do
	  end do
 
	  ! Intramolecular distribution
	  do n1=-pgrid,pgrid
	    do n2=-pgrid,pgrid
	      do n3=-pgrid,pgrid
		!pdensintra(sp1,n1,n2,n3)=pdensintra(sp1,n1,n2,n3)/nframes/s_nmols(sp1)/(pdelta**3)
		pdensintra(sp1,n1,n2,n3)=pdensintra(sp1,n1,n2,n3)/norm2/(pdelta**3)
	      end do
	    end do
	  end do

	  ! Average molecule
	  do n=1,s_natoms(sp1)
	    do m=1,3
	      avggeom(sp1,n,m)=avggeom(sp1,n,m)/norm2
	    end do
	  end do

	end do

	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "Problem with HISTORY file - failed to read header."
	goto 999
800	write(0,*) "End of unformatted HISTORY file found."
	write(15,"(A)") "End of unformatted HISTORY file found."
801	write(0,"(A,I1)") " Central species = ",centresp
	do sp1=1,nspecies
	  write(0,"(A,I1,A,A)") " Species ",sp1,", ",s_name(sp1)
	  write(0,"(A12,I9,A,I9,A)") "Expected : ",spexp(sp1)," over ",nframes," frames."
	  write(0,"(A12,I9)") "Found : ",pdensn(sp1)
	  write(0,"(A12,I9,A)") "Caught : ",pdensa(sp1)," (in grid)"
	  if (molmap.and.(sp1.eq.centresp)) write(0,"(A12,I9,A)") "Selected : ",molcount," (using mapfile)"
	  write(15,"(A,I1,A,A)") " Species ",sp1,", ",s_name(sp1)
	  write(15,"(A12,I9,A,I9,A)") "Expected : ",spexp(sp1)," over ",nframes," frames."
	  write(15,"(A12,I9)") "Found : ",pdensn(sp1)
	  write(15,"(A12,I9,A)") "Caught : ",pdensa(sp1)," (in grid)"
	  if (molmap.and.(sp1.eq.centresp)) write(15,"(A12,I9,A)") "Selected : ",molcount," (using mapfile)"
	end do
	write(0,"(3(A,I3),A,F4.1,A)") "Grid = ",(grid*2+1),"x",(grid*2+1),"x", &
	  & (grid*2+1)," ( x ",delta," Angstrom )"
	write(15,"(3(A,I3),A,F4.1,A)") "Grid = ",(grid*2+1),"x",(grid*2+1),"x", &
	  & (grid*2+1)," ( x ",delta," Angstrom )"
	write(0,"(\,A,\)") "-- Writing output files..."

	do sp1=1,nspecies
	  ! -- Probability density of this molecule about the central species
	  resfile=basename(1:baselen)//"pdens"//CHAR(48+centresp)//CHAR(48+sp1)
	  open(unit=9,file=resfile,form="formatted")
	  do n=-grid,grid
	    do m=-grid,grid
	      do o=-grid,grid
	  	write(9,"(F12.8)") pdens(sp1,n,m,o)
	      end do
	    end do
	  end do
	  close(9)

	  ! -- Average molecule
	  resfile=basename(1:baselen)//"avgmol"//char(48+sp1)
	  open(unit=9,file=resfile,form="formatted")
	  ! Write out pdb files of the average molecules
900	  FORMAT ('ATOM',4x,I3,2x,A2,15x,3(F7.3,2x))
901	  FORMAT ('ATOM',4x,I3,2x,A2,15x,3(F12.5,2x))
	  do n=1,s_natoms(sp1)
	    write(9,900) n,atmname(s_start(sp1)-1+n)(1:2),(avggeom(sp1,n,m),m=1,3)
	  end do
	  close(9)

	  ! -- Intramolecular probability distribution
	  resfile=basename(1:baselen)//"intra"//char(48+sp1)
	  open(unit=9,file=resfile,form="formatted")
902	  FORMAT (3F7.4,F12.8)
	  do n=-pgrid,pgrid
	    do m=-pgrid,pgrid
	      do o=-pgrid,pgrid
	        !write(9,*) n*pdelta,m*pdelta,o*pdelta,pdensintra(sp1,n,m,o)
	        write(9,"(F12.8)") pdensintra(sp1,n,m,o)
	      end do
	    end do
	  end do
	  close(9)

	end do

	write(0,*) "Finished!"
	write(15,"(A)") "Finished!"
999	close(10)
	close(13)

	end program pdens3d
