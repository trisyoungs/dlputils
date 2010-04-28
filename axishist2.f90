!	** axisdens2.f90 **
!	Calculate the orientation of other species around a central one, binned into angular regions

	program axishist2
	use dlprw; use utility
	implicit none
	real*8, parameter :: radcon = 57.29577951d0
	integer, allocatable :: ahist(:,:,:,:,:), acc(:,:)
	integer :: baselen,nargs,nframesused
	character*80 :: hisfile,outfile,basename,resfile,temp,flagfile,altheaderfile
	character*8 :: discard
	logical :: altheader = .FALSE.
	integer :: success,binx,biny,binz,nanglebins,anglexbin,angleybin
	integer :: iatm, n, nframes, k, l, m, o, p, i
	real*8 :: px,py,pz, ax, ay, az, anglebin
	real*8 :: mindist, maxdist, dist, tx,ty,tz, total
	integer :: sp1,m1,sp2,m2,startf,endf
	integer :: iargc

	write(0,*) "*** axishist"

	nargs = iargc()
	if (nargs.lt.4) then
	  write(*,"(a)") "Usage: axishist2 <HIStory file> <OUTput file> <sp1> <sp2>"
	  write(*,"(a)") "        [-axis sp x1 x2 y1 y2]  Atoms to use for axis calculation in species sp"
	  write(*,"(a)") "        [-grid npoints]         Grid points in each direction for prob. densities (default = 20)"
	  write(*,"(a)") "        [-delta spacing]        Spacing between grid points (default = 0.5 Angstrom)"
	  write(*,"(a)") "        [-start frameno]        Trajectory frame to start calculations (default = 1)"
	  write(*,"(a)") "        [-end frameno]          Trajectory frame to end calculations on (default = last)"
	  write(*,"(a)") "        [-mindist r]            Minimum separation allowed between molecules (default = 0.0)"
	  write(*,"(a)") "        [-maxdist r]            Maximum separation allowed between molecules (default = 10000.0)"
	  write(*,"(a)") "        [-anglebin d]           Size of angle bin to use (default = 15deg)"
	  write(*,"(a)") "        [-header file]          Use specified file to get header"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
	  call getarg(3,temp); read(temp,"(I4)") sp1
	  call getarg(4,temp); read(temp,"(I4)") sp2
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
	  basename="adensresults."
	  baselen=13
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=15,file=basename(1:baselen)//"out",form="formatted",status="replace")

	! Set some variable defaults before we read in any command-line arguments
	call alloc_axis
	mindist = 0.0
	maxdist = 10000.0
	anglebin = 15.0
	startf = 1
	endf = 0

	n = 4
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-axis") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") m
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(m,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(m,2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(m,3)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(m,4)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",m," calculated from: X=",aa(m,1),"->", &
	        & aa(m,2),", Y=0.5X->0.5(r(",aa(m,3),")->r(",aa(m,4),"))"
	      axisdefined(m) = .true.
	    case ("-start")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") startf
	      write(0,"(A,I5)") "Starting frame = ",startf
	    case ("-end")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") endf
	      write(0,"(A,I5)") "End frame = ",endf
	    case ("-mindist")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") mindist
	      write(0,"(A,f8.4)") "Minimum separation between molecules = ",mindist
	    case ("-maxdist")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") maxdist
	      write(0,"(A,f8.4)") "Maximum separation between molecules = ",maxdist
	    case ("-anglebin")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") anglebin
	      write(0,"(A,f8.4)") "Size of angular bin = ",anglebin
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
	do n=1,nspecies
	  ! Check that molecules have had their axes defined
	  if (axisdefined(n)) then
	    write(15,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Axis for species ",n," calculated from : X=",aa(n,1),"->", &
	      & aa(n,2),", Y=0.5X->0.5(",aa(n,3),"->",aa(n,4),")"
	  else
	    if (s_natoms(n).gt.2) stop "Axes must be defined for all molecular species."
	  end if
	end do
	write(15,"(A,I4)") "Central species = ",sp1
	write(15,"(A,I4)") "Surrounding species = ",sp2
	write(15,"(A,I5,A,I5,A)") "Frame range = ",startf," to ",endf," (0=last)"

	! Allocate and clear the arrays before we start...
	write(0,*) "Allocating and clearing initial arrays..."
	nanglebins = int(180.0 / anglebin) + 1
	allocate(ahist(nanglebins,nanglebins,3,3,-100:100))
	allocate(acc(nanglebins,nanglebins))
	ahist = 0
	acc = 0

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
100	nframes = 0
	nframesused = 0
102	success=readframe()
	IF (success.eq.1) goto 120  ! End of file encountered....
	IF (success.eq.-1) goto 119  ! File error....
	nframes=nframes+1
	if (MOD(nframes,100).eq.0) write(0,"(i6)") nframes
	if (nframes.LT.startf) goto 102

	nframesused = nframesused + 1

	! Generate all molecular axes
	call genaxis

	! Calculate histogram
	do m1=1,s_nmols(sp1)      ! Loop over all molecules of species 1...

	  do m2=1,s_nmols(sp2)

	    ! Don't consider same molecule with itself
	    if ((m1.eq.m2).and.(sp1.eq.sp2)) cycle

	    ! Get minimum image position of second centre with first
	    call pbc(axisox(sp2,m2),axisoy(sp2,m2),axisoz(sp2,m2),axisox(sp1,m1),axisoy(sp1,m1),axisoz(sp1,m1),tx,ty,tz)
	    tx = tx - axisox(sp1,m1)
	    ty = ty - axisoy(sp1,m1)
	    tz = tz - axisoz(sp1,m1)
	    dist = sqrt(tx*tx + ty*ty + tz*tz)
	    if (dist.lt.mindist) cycle
	    if (dist.gt.maxdist) cycle

	    ! Determine angular bins
	    ! Vectors COM-X an COM-Y are the cartesian unit vectors (1,0,0) and (0,1,0)
	    ! tx,ty,tz contain the distance vector between the axis origins of sp1/m1 and sp2/m2
	    tx = tx / dist
	    ty = ty / dist
	    tz = tz / dist
	    anglexbin = abs( acos( 1.0*tx )*radcon / anglebin + 1 )
	    angleybin = abs( acos( 1.0*ty )*radcon / anglebin + 1 )
	!write(0,*) "anglex = ", acos( 1.0*tx )*radcon
	!write(0,*) "angley = ", acos( 1.0*ty )*radcon
	!write(0,*) anglexbin, angleybin
	    

	    ! This molecule is within the limits defined, so rotate its axes into the reference frame one by one, binning along the way
            ax = axisx(sp2,m2,1)*axisx(sp1,m1,1) + axisx(sp2,m2,2)*axisx(sp1,m1,2) + axisx(sp2,m2,3)*axisx(sp1,m1,3)
            ay = axisx(sp2,m2,1)*axisy(sp1,m1,1) + axisx(sp2,m2,2)*axisy(sp1,m1,2) + axisx(sp2,m2,3)*axisy(sp1,m1,3)
            az = axisx(sp2,m2,1)*axisz(sp1,m1,1) + axisx(sp2,m2,2)*axisz(sp1,m1,2) + axisx(sp2,m2,3)*axisz(sp1,m1,3)
	    binx = ax / 0.01
	    biny = ay / 0.01
	    binz = az / 0.01
	    ahist(anglexbin,angleybin,1,1,binx) = ahist(anglexbin,angleybin,1,1,binx) + 1
	    ahist(anglexbin,angleybin,1,2,biny) = ahist(anglexbin,angleybin,1,2,biny) + 1
	    ahist(anglexbin,angleybin,1,3,binz) = ahist(anglexbin,angleybin,1,3,binz) + 1

            ax = axisy(sp2,m2,1)*axisx(sp1,m1,1) + axisy(sp2,m2,2)*axisx(sp1,m1,2) + axisy(sp2,m2,3)*axisx(sp1,m1,3)
            ay = axisy(sp2,m2,1)*axisy(sp1,m1,1) + axisy(sp2,m2,2)*axisy(sp1,m1,2) + axisy(sp2,m2,3)*axisy(sp1,m1,3)
            az = axisy(sp2,m2,1)*axisz(sp1,m1,1) + axisy(sp2,m2,2)*axisz(sp1,m1,2) + axisy(sp2,m2,3)*axisz(sp1,m1,3)
	    binx = ax / 0.01
	    biny = ay / 0.01
	    binz = az / 0.01
	    ahist(anglexbin,angleybin,2,1,binx) = ahist(anglexbin,angleybin,2,1,binx) + 1
	    ahist(anglexbin,angleybin,2,2,biny) = ahist(anglexbin,angleybin,2,2,biny) + 1
	    ahist(anglexbin,angleybin,2,3,binz) = ahist(anglexbin,angleybin,2,3,binz) + 1

            ax = axisz(sp2,m2,1)*axisx(sp1,m1,1) + axisz(sp2,m2,2)*axisx(sp1,m1,2) + axisz(sp2,m2,3)*axisx(sp1,m1,3)
            ay = axisz(sp2,m2,1)*axisy(sp1,m1,1) + axisz(sp2,m2,2)*axisy(sp1,m1,2) + axisz(sp2,m2,3)*axisy(sp1,m1,3)
            az = axisz(sp2,m2,1)*axisz(sp1,m1,1) + axisz(sp2,m2,2)*axisz(sp1,m1,2) + axisz(sp2,m2,3)*axisz(sp1,m1,3)
	    binx = ax / 0.01
	    biny = ay / 0.01
	    binz = az / 0.01
	    ahist(anglexbin,angleybin,3,1,binx) = ahist(anglexbin,angleybin,3,1,binx) + 1
	    ahist(anglexbin,angleybin,3,2,biny) = ahist(anglexbin,angleybin,3,2,biny) + 1
	    ahist(anglexbin,angleybin,3,3,binz) = ahist(anglexbin,angleybin,3,3,binz) + 1

	    acc(anglexbin,angleybin) = acc(anglexbin,angleybin) + 1

	  end do

	end do

	! Next frame (or finish)
116	if (nframes.EQ.endf) goto 120
	goto 102

119	write(0,*) "HISTORY file ended prematurely!"
	write(0,"(A,I5,A)") "Managed ",nframesused," frames before error."
	write(15,"(A)") "HISTORY file ended prematurely!"
	write(15,"(A,I5,A)") "Managed ",nframesused," frames before error."
120	write(0,*) "Finished."
	write(15,"(A)") "Finished."
	write(15,"(A,I5,A)") "Averages will be taken over ",nframesused," frames."

	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "Problem with HISTORY file - failed to read header."
	goto 999
800	write(0,*) "End of unformatted HISTORY file found."
	write(15,"(A)") "End of unformatted HISTORY file found."
801	write(0,"(A)") "-- Writing output files..."

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).EQ.".") THEN
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) THEN
	  basename="rdfresults."
	  baselen=11
	ELSE
	  basename=hisfile(1:baselen)
	endif

	! Write accumulation data
        resfile=basename(1:baselen)//"ahist2.acc"
	open(unit=9, file=resfile, form="formatted", status="replace")
        resfile=basename(1:baselen)//"ahist2.acc.gp"
	open(unit=10, file=resfile, form="formatted", status="replace")
	write(9,"(6x,30(f5.1,1x))") (m*anglebin, m=0,nanglebins-1)
	total = sum(acc)
	do n=1,nanglebins+1
	  write(9,"(f5.1,1x,30(f5.3,1x))") (n-1)*anglebin, (acc(n,m)/total,m=1,nanglebins+1)
	  if (n.le.nanglebins) then
	    do m=1,nanglebins
	      write(10,*) acc(n,m)/total
	    end do
	    write(10,*) ""
	  end if
	end do
	close(9)
	close(10)

	! Write unnormalised histograms
	total = nframesused*1.0
	do n=1,nanglebins
	  do m=1,nanglebins
	    resfile=basename(1:baselen)//"ahist"//CHAR(48+sp1)//"_"//CHAR(48+sp2)//".x"//CHAR(48+(n/10)) &
		& //CHAR(48+MOD(n,10))//"_y"//CHAR(48+(m/10))//CHAR(48+MOD(m,10))
	    open(unit=9, file=resfile, form="formatted", status="replace")
	    write(9,"(10a10)") "#D","XX","XY","XZ","YX","YY","YZ","ZX","ZY","ZZ"
	    ! Correct central bin
	    do o=1,3
	      do p=1,3
	        ahist(n,m,o,p,0) = ahist(n,m,o,p,0) / 2
	      end do
	    end do
	    ! Write data
	    do i=-100,100
	      write(9,"(10f10.4)") real(i)*0.01,(ahist(n,m,1,p,i)/total,p=1,3),(ahist(n,m,2,p,i)/total,p=1,3),(ahist(n,m,3,p,i)/total,p=1,3)
	    end do
	    close(9)
	  end do
	end do

	write(0,*) "Finished!"
	write(15,"(A)") "Finished!"
999	close(10)
	close(11)

	end program axishist2
