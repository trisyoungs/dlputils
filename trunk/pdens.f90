!	** pdens.f90 **
!	Program to calculate the 3-dimensional distribution of different species about
!	one species' centres.
!	Updated version - rotates all molecules about target molecule first, then performs minimum
!	image before binning

	program pdens2
	use dlprw; use utility
	implicit none
	real*8, allocatable :: pdens(:,:,:,:)
	real*8, allocatable :: avggeom(:,:,:)		! Average species coordinates
	real*8, allocatable :: pdensintra(:,:,:,:)	! Intramolecular species distributions
	real*8, allocatable :: rrot(:,:)		! Rotated coordinates of axis centres
	integer, allocatable :: nfound(:), ncaught(:)	! Total and 'binned' mols
	integer, allocatable :: molflags(:)		! Per-frame map of sp1 mols to include in averaging
	integer, allocatable :: spexp(:)		! Expected species numbers
	integer :: grid = 20, pgrid = 50		! Grid in each direction for 3ddist and pdens
	real*8 :: delta = 0.5, pdelta = 0.15		! Default grid spacings
	logical :: molmap = .FALSE.			! Whether we're using a file of mol flags
	logical :: nopdens = .FALSE., nointra = .FALSE.	! Flags to prohibit various calculations
	integer :: nmapframes, baselen,nargs,nframesused
	character*80 :: hisfile,outfile,basename,resfile,temp,flagfile,altheaderfile
	logical :: altheader = .FALSE.
	integer :: success,n1,n2,n3
	integer :: n, nframes, m, o, p
	real*8 :: px, py, pz, mindist, maxdist, dist, tx, ty, tz
	integer :: centresp,sp1,m1,sp2,m2,startf,endf,n3dcentres
	integer :: iargc

	write(0,*) "*** pdens"

	nargs = iargc()
	if (nargs.LT.2) then
	  write(*,"(a)") "Usage: pdens <HIStory file> <OUTput file> [-options]"
	  write(*,"(a)") "        [-centre sp]            Set species sp to be the central one"
	  write(*,"(a)") "        [-axis sp x1 x2 y1 y2]  Atoms to use for axis calculation in species sp"
	  write(*,"(a)") "        [-grid npoints]         Grid points in each direction for prob. densities (default = 20)"
	  write(*,"(a)") "        [-delta spacing]        Spacing between grid points (default = 0.5 Angstrom)"
	  write(*,"(a)") "        [-pgrid npoints]        Grid points in each direction for species densities (default = 50)"
	  write(*,"(a)") "        [-pdelta spacing]       Spacing between pgrid points (default = 0.15 Angstrom)"
	  write(*,"(a)") "        [-start frameno]        Trajectory frame to start calculations (default = 1)"
	  write(*,"(a)") "        [-end frameno]          Trajectory frame to end calculations on (default = last)"
	  write(*,"(a)") "        [-molmap <file>]        Formatted file (I5,1X,n(I2,1x) of frame,mols to use of sp1"
	  write(*,"(a)") "        [-nopdens,-nointra]     Prohibit parts of the calculation"
	  write(*,"(a)") "        [-mindist r]            Minimum separation allowed between molecules (default = 0.0)"
	  write(*,"(a)") "        [-maxdist r]            Maximum separation allowed between molecules (default = 10000.0)"
	  write(*,"(a)") "        [-header file]          Use specified file to get header"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
	end if
	n = 10
	write(0,"(A,A)") "History file : ",hisfile
	!call openhis(hisfile,n)
	write(0,"(A,A)") " Output file : ",outfile
	if (outinfo(outfile,1).eq.-1) goto 798

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  IF (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 50
	  endIF
	end do
50	if (baselen.eq.-1) then
	  basename="3ddistresults."
	  baselen=14
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=15,file=basename(1:baselen)//"out",form="formatted",status="replace")

	! Set some variable defaults before we read in any command-line arguments
	call alloc_axis()
	mindist = 0.0
	maxdist = 10000.0
	centresp = 1
	startf = 1
	endf = 0
	n3dcentres = 0

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
	      axisdefined(sp1) = .true.
	      if (aa(sp1,1).eq.aa(sp1,2)) then
		write(0,*) " --- X-axis atom ids are identical - atom will be used as centre but NO axes will be defined."
	        axisdefined(sp1) = .false.
		axisusecom(sp1) = .false.
	      end if
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
	      write(0,"(A,I5)") "Starting frame = ",startf
	    case ("-end")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") endf
	      write(0,"(A,I5)") "End frame = ",endf
	    case ("-molmap")
	      n = n + 1; call getarg(n,flagfile)
	      write(0,"(A,A)") "Using flags for species 1 molecules from : ",flagfile
	      open(unit=20,file=flagfile,form="formatted",status="old")
	      molmap = .TRUE.
	    case ("-nopdens")
	      nopdens = .TRUE.
	      write(0,"(A,I4)") "Intermolecular probability distributions will not be computed."
	      write(15,"(A,I4)") "Intermolecular probability distributions will not be computed."
	    case ("-nointra")
	      nointra = .TRUE.
	      write(0,"(A,I4)") "Intramolecular probability distributions will not be computed."
	      write(15,"(A,I4)") "Intramolecular probability distributions will not be computed."
	    case ("-mindist")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") mindist
	      write(0,"(A,f8.4)") "Minimum separation between molecules = ",mindist
	      write(15,"(A,f8.4)") "Minimum separation between molecules = ",mindist
	    case ("-maxdist")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") maxdist
	      write(0,"(A,f8.4)") "Maximum separation between molecules = ",maxdist
	      write(15,"(A,f8.4)") "Maximum separation between molecules = ",maxdist
	    case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
              altheader = .TRUE.
	    case default
	      write(0,*) "Unrecognised command line option : ",temp
	      stop
	  end select
	end do
	
	! Print out a summary of the control variables to the output file.
	write(15,"(A,A)") "Input file: ",hisfile
	write(15,"(A,I3)") "Molecular species in file : ",nspecies
	do sp1=1,nspecies
	  ! Check that molecules have had their axes defined
	  if (axisdefined(sp1)) then
	    write(15,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Axis for species ",sp1," calculated from : X=",aa(sp1,1),"->", &
	      & aa(sp1,2),", Y=0.5X->0.5(",aa(sp1,3),"->",aa(sp1,4),")"
	  else if (sp1.eq.centresp) then
	    stop "A proper set of axes must be defined for the central species."
	  end if
	end do
	write(15,"(A,I4)") "Grid (intermolecular) points in each XYZ = ",grid
	write(15,"(A,I4)") "PGrid (intramolecular) points in each XYZ = ",pgrid
	write(15,"(A,F6.3)") "Grid spacing = ",delta
	write(15,"(A,F6.3)") "PGrid spacing = ",pdelta
	write(15,"(A,I4)") "Central species = ",centresp
	write(15,"(A,I5,A,I5,A)") "Frame range = ",startf," to ",endf," (0=last)"

	! Probability density arrays
	allocate(pdens(nspecies,-grid:grid,-grid:grid,-grid:grid))
	allocate(ncaught(nspecies))
	allocate(nfound(nspecies))
	allocate(pdensintra(nspecies,-pgrid:pgrid,-pgrid:pgrid,-pgrid:pgrid))
	allocate(spexp(nspecies))
	! Rotated coordinates
	allocate(rrot(maxmols,3))
	! Average species etc.
	allocate(avggeom(nspecies,maxatoms,3))
	if (molmap) allocate(molflags(s_nmols(centresp)))

	! Clear the arrays before we start...
	write(0,*) "Clearing initial arrays..."
	pdens = 0.0
	ncaught = 0
	nfound = 0
	pdensintra = 0.0
	spexp = 0
	avggeom = 0.0

	! Read the header of the history file...
        call openhis(hisfile,10)
        if (readheader().EQ.-1) then
          if (altheader) then
            write(0,*) "Restarted trajectory:"
            close(dlpun_his)
            call openhis(altheaderfile,10)
            if (readheader().EQ.-1) goto 799
            close(dlpun_his)
            call openhis(hisfile,10)
          else
            goto 799
          end if
        end if

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
	nframes = 0
	nframesused = 0
	nmapframes = 0
	! If we're using a mapfile, read it here
101	if (molmap) then
	  read(20,"(I5,1x,40(I2,1x))",end=118,err=117) m1,(molflags(m2),m2=1,s_nmols(centresp))
	  ! if (m1.NE.nframes) stop "Mapfile missed a frame number!"
	end if
102	success=readframe()
	IF (success.eq.1) goto 120  ! End of file encountered....
	IF (success.lt.0) goto 119  ! File error....
	nframes=nframes+1

	if (MOD(nframes,100).eq.0) then
	  if (molmap) then
	    write(0,"(I6,2x,'(',I6,',',I10,')')") nframes,nmapframes,n3dcentres
	  else
	    write(0,"(i6,2x,'(',I10,')')") nframes,n3dcentres
	  end if
	end if
	if (nframes.LT.startf) goto 102

	if ((molmap).and.(m1.ne.nframes)) goto 102
	if (molmap) nmapframes = nmapframes + 1

	nframesused = nframesused + 1

	! Generate all molecular axes
	call genaxis

	!
	! Calculate 3D distributions
	!
	if (nopdens) goto 115
	sp1 = centresp
	do m1=1,s_nmols(sp1)      ! Loop over all molecules of species 1...
	  ! If we're using a mapfile, decide whether to include this molecule
	  if (molmap) then
	    if (molflags(m1).eq.0) cycle 
	  end if
	  n3dcentres = n3dcentres + 1	! Normalisation counter
	  do sp2=1,nspecies
	    do m2=1,s_nmols(sp2)

	      ! Calculate minimum image position of molecule m2 about 0,0,0
	      call pbc(axisox(sp2,m2),axisoy(sp2,m2),axisoz(sp2,m2),axisox(sp1,m1),axisoy(sp1,m1),axisoz(sp1,m1),tx,ty,tz)
	      px = tx - axisox(sp1,m1)
	      py = ty - axisoy(sp1,m1)
	      pz = tz - axisoz(sp1,m1)
	      tx = px*axisx(sp1,m1,1) + py*axisx(sp1,m1,2) + pz*axisx(sp1,m1,3)
	      ty = px*axisy(sp1,m1,1) + py*axisy(sp1,m1,2) + pz*axisy(sp1,m1,3)
	      tz = px*axisz(sp1,m1,1) + py*axisz(sp1,m1,2) + pz*axisz(sp1,m1,3)

	      ! Get relative coordinates (i.e. centre new coords about 0,0,0)
	      !tx = axisox(sp2,m2) - axisox(sp1,m1)
	      !ty = axisoy(sp2,m2) - axisoy(sp1,m1)
	      !tz = axisoz(sp2,m2) - axisoz(sp1,m1)

	      ! Rotate relative coordinates into local frame of molecule m1
	      !px = tx*axisx(sp1,m1,1) + ty*axisx(sp1,m1,2) + tz*axisx(sp1,m1,3)
	      !py = tx*axisy(sp1,m1,1) + ty*axisy(sp1,m1,2) + tz*axisy(sp1,m1,3)
	      !pz = tx*axisz(sp1,m1,1) + ty*axisz(sp1,m1,2) + tz*axisz(sp1,m1,3)

	      ! Calculate minimum image position of molecule m2 about 0,0,0
	      !tx = px - cell(1)*NINT(px/cell(1))
	      !ty = py - cell(5)*NINT(py/cell(5))
	      !tz = pz - cell(9)*NINT(pz/cell(9))

	      ! Check minimum / maximum separation
	      dist = sqrt(tx*tx + ty*ty + tz*tz)
	      if (dist.lt.mindist) cycle
	      if (dist.gt.maxdist) cycle

	      ! Calculate integer position in grid
	      n1=NINT(tx/delta)
	      n2=NINT(ty/delta)
	      n3=NINT(tz/delta)

	      ! If any of the n's are over grid, then only increase the counter
	      if (MAX(ABS(n1),MAX(ABS(n2),ABS(n3))).GT.grid) then
	        nfound(sp2) = nfound(sp2) + 1      ! No position, just count it....
	      else
		if ((sp1.eq.sp2).and.(m1.eq.m2)) then
		  ! Same molecule, so don't add at all
		else
	          ! Check to make sure the same molecule isn't being consider with itself  XXXXX
	          pdens(sp2,n1,n2,n3) = pdens(sp2,n1,n2,n3) + 1
	          nfound(sp2) = nfound(sp2) + 1
	          ncaught(sp2) = ncaught(sp2) + 1
		end if
	      end if
	      p=p+s_natoms(sp2)
	    end do
	  end do
	end do

	! Calculate average molecules and intramolecular distributions.
115	if (nointra) goto 116
	do sp1=1,nspecies
	  p=s_start(sp1)
	  do m1=1,s_nmols(sp1)
	    ! If we're using a mapfile, decide whether to include this molecule
	    if ((sp1.eq.centresp).AND.(molmap)) then
	      if (molflags(m1).eq.0) then
		p = p+s_natoms(sp1)
		cycle
	      end if
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
	      if (MAX(ABS(n1),MAX(ABS(n2),ABS(n3))).lt.pgrid) then
	        pdensintra(sp1,n1,n2,n3)=pdensintra(sp1,n1,n2,n3)+1
	      else
		write(0,*) "Large distance : n1,n2,n3 = ",n1,n2,n3
	      end if
	    end do
	    p=p+s_natoms(sp1)
	  end do
	end do

	! Next frame (or finish)
116	if (nframes.EQ.endf) goto 120
	goto 101

117	write(0,*) "Error reading mapfile."
	write(0,*) "Selected ",nmapframes," from trajectory before error."
	goto 120
118	write(0,*) "Reached end of mapfile."
	write(0,*) "Selected ",nmapframes," from trajectory."
	goto 120
119	write(0,*) "HISTORY file ended prematurely!"
	write(0,"(A,I5,A)") "Managed ",nframesused," frames before error."
	write(15,"(A)") "HISTORY file ended prematurely!"
	write(15,"(A,I5,A)") "Managed ",nframesused," frames before error."
120	write(0,*) "Finished."
	write(15,"(A)") "Finished."
	write(15,"(A,I5,A)") "Averages will be taken over ",nframesused," frames."

	! ### Normalise the data

	write(0,"(a,i1,a,i7,a,i7,a)") "Central species consisted ",n3dcentres," molecules"

	! Set normalisation over frames and central species molecules
	do sp1=1,nspecies

	  ! Calculate expected species numbers
	  spexp(sp1) = s_nmols(centresp) * s_nmols(sp1) * nframesused
	  if (centresp.eq.sp1) spexp(sp1) = spexp(sp1) - s_nmols(sp1) * nframesused

	  ! Species density about central species
	  do n1=-grid,grid
	    do n2=-grid,grid
	      do n3=-grid,grid
	        ! pdens(sp1,n1,n2,n3)=pdens(sp1,n1,n2,n3)/norm1/(delta**3)
	        pdens(sp1,n1,n2,n3)=pdens(sp1,n1,n2,n3)/n3dcentres/(delta**3)
	      end do
	    end do
	  end do
 
	  ! Intramolecular distribution
	  do n1=-pgrid,pgrid
	    do n2=-pgrid,pgrid
	      do n3=-pgrid,pgrid
		!pdensintra(sp1,n1,n2,n3)=pdensintra(sp1,n1,n2,n3)/nframes/s_nmols(sp1)/(pdelta**3)
		if (sp1.eq.centresp) then
		  pdensintra(sp1,n1,n2,n3)=pdensintra(sp1,n1,n2,n3)/n3dcentres/(pdelta**3)
		else
		  pdensintra(sp1,n1,n2,n3)=pdensintra(sp1,n1,n2,n3)/(nframesused*s_nmols(sp1))/(pdelta**3)
		end if
	      end do
	    end do
	  end do

	  ! Average molecule
	  do n=1,s_natoms(sp1)
	    do m=1,3
	      if (sp1.eq.centresp) then
	        avggeom(sp1,n,m)=avggeom(sp1,n,m)/n3dcentres
	      else
	        avggeom(sp1,n,m)=avggeom(sp1,n,m)/(nframesused*s_nmols(sp1))
	      end if
	    end do
	  end do

	end do

	goto 801

	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "Problem with HISTORY file - failed to read header."
	goto 999
	write(0,*) "End of unformatted HISTORY file found."
	write(15,"(A)") "End of unformatted HISTORY file found."
801	write(0,"(A,I1)") " Central species = ",centresp
	do sp1=1,nspecies
	  write(0,"(A,I1,A,A)") " Species ",sp1,", ",s_name(sp1)
	  write(0,"(A12,I9,A,I9,A)") "Expected : ",spexp(sp1)," over ",nframesused," frames."
	  write(0,"(A12,I9)") "Found : ",nfound(sp1)
	  write(0,"(A12,I9,A)") "Caught : ",ncaught(sp1)," (in grid)"
	  if (sp1.eq.centresp) write(0,"(A12,I9)") "Selected : ",n3dcentres
	  write(15,"(A,I1,A,A)") " Species ",sp1,", ",s_name(sp1)
	  write(15,"(A12,I9,A,I9,A)") "Expected : ",spexp(sp1)," over ",nframesused," frames."
	  write(15,"(A12,I9)") "Found : ",nfound(sp1)
	  write(15,"(A12,I9,A)") "Caught : ",ncaught(sp1)," (in grid)"
	  if (sp1.eq.centresp) write(15,"(A12,I9)") "Selected : ",n3dcentres
	end do
	write(0,"(3(A,I3),A,F4.1,A)") "Grid = ",(grid*2+1),"x",(grid*2+1),"x", &
	  & (grid*2+1)," ( x ",delta," Angstrom )"
	write(15,"(3(A,I3),A,F4.1,A)") "Grid = ",(grid*2+1),"x",(grid*2+1),"x", &
	  & (grid*2+1)," ( x ",delta," Angstrom )"
	write(0,"(A)") "-- Writing output files..."

	do sp1=1,nspecies

	  if (.not.nopdens) then
	    ! -- Probability density of this molecule about the central species
	    resfile=basename(1:baselen)//CHAR(48+centresp)//CHAR(48+sp1)//".pdens"
	    open(unit=9,file=resfile,form="formatted")
	    write(9,*) 2*grid+1, 2*grid+1, 2*grid+1
	    write(9,"(9f8.4)") delta,0.0,0.0,0.0,delta,0.0,0.0,0.0,delta
	    write(9,"(3f10.4)") -grid*delta,-grid*delta,-grid*delta
	    write(9,*) "zyx"
	    do n=-grid,grid
	      do m=-grid,grid
	        do o=-grid,grid
	    	write(9,"(f12.8)") pdens(sp1,n,m,o)
	        end do
	      end do
	    end do
	    close(9)
	  end if

	  if (.not.nointra) then
	    ! -- Average molecule
	    resfile=basename(1:baselen)//"avg"//char(48+sp1)//".xyz"
	    open(unit=9,file=resfile,form="formatted")
	    ! Write out pdb files of the average molecules
900	    FORMAT (a2,5x,3(F9.5,2x))
	    write(9,*) s_natoms(sp1)
	    write(9,*) "Average: ", s_name(sp1)
	    do n=1,s_natoms(sp1)
	      write(9,900) atmname(s_start(sp1)-1+n)(1:2),(avggeom(sp1,n,m),m=1,3)
	    end do
	    close(9)

	    ! -- Intramolecular probability distribution
	    resfile=basename(1:baselen)//"intra"//char(48+sp1)//".pdens"
	    open(unit=9,file=resfile,form="formatted")
	    write(9,*) 2*pgrid+1, 2*pgrid+1, 2*pgrid+1
	    write(9,"(9f8.4)") pdelta,0.0,0.0,0.0,pdelta,0.0,0.0,0.0,pdelta
	    write(9,"(3f10.4)") -pgrid*pdelta,-pgrid*pdelta,-pgrid*pdelta
	    write(9,*) "zyx"
	    do n=-pgrid,pgrid
	      do m=-pgrid,pgrid
	        do o=-pgrid,pgrid
	          !write(9,*) n*pdelta,m*pdelta,o*pdelta,pdensintra(sp1,n,m,o)
	          write(9,"(F12.8)") pdensintra(sp1,n,m,o)
	        end do
	      end do
	    end do
	    close(9)
	  end if

	end do

	write(0,*) "Finished!"
	write(15,"(A)") "Finished!"
999	close(10)
	close(13)

	end program pdens2
