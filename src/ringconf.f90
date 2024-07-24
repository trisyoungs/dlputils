	! ringconf
	! Calculates the spherical polar coordinate description of 6-rings, using the method
	! of Berces, Whitfield and Nukuda, described in Tetrahedron, 57 (3), 477 (2001)
	! Many thanks to them for providing the original Fortran codes on which parts of this code are based

	module conformerData
	  ! Conformer Data
	  integer, parameter :: CHAIR=1, BOAT=2, TWISTBOAT=3
	  real*8, parameter :: radcon = 57.29577951d0
	  real*8, parameter :: pi = 3.14159265358979d0

	  character*3 :: basicnames(38)
	  integer :: ideal_torsions(3,6), redundant_torsions(3,6), basic_torsions(38,6)
	  real*8 :: ideal_norms(3), redundant_norms(3), basic_norms(38)
	contains

	subroutine setupConformerData()
	implicit none

	integer :: n
	real*8 :: tempv(6)

	! Set ideal, redundancy, and basic torsion vectors.
	! Redundant torsions R4, R5, and R6 map to redundant_torsions(1-3,:)
	ideal_torsions(CHAIR,:)     = (/ 60, -60, 60, -60, 60, -60 /)	! 1C4
	ideal_torsions(BOAT,:)      = (/ 0, 60, -60, 0, 60, -60 /)	! 14B
	ideal_torsions(TWISTBOAT,:) = (/ 60, -30, -30, 60, -30, -30 /)	! OS2
	redundant_torsions(TWISTBOAT,:) = (/ 60, 30, -30, -60, -30, 30 /)	! Twist
	redundant_torsions(BOAT,:)     = (/ 0, 60, 60, 0, -60, -60 /)	! Boat
	redundant_torsions(CHAIR,:)     = (/ 60, 60, 60, 60, 60, 60 /)	! Chair

	basic_torsions(1,:) = (/  60, -60,  60, -60,  60, -60 /)	! 1C4
	basic_torsions(2,:) = (/ -60,  60, -60,  60, -60,  60 /)	! 4C1

	basic_torsions(3,:)  = (/   0,  60, -60,   0,  60, -60 /)	! "1,4B"
	basic_torsions(4,:)  = (/  60,   0, -60,  60,   0, -60 /)	! "B2,5"
	basic_torsions(5,:)  = (/  60, -60,   0,  60, -60,   0 /)	! "O,3B"
	basic_torsions(6,:)  = (/   0, -60,  60,   0, -60,  60 /)	! "B1,4"
	basic_torsions(7,:)  = (/ -60,   0,  60, -60,   0,  60 /)	! "2,5B"
	basic_torsions(8,:)  = (/ -60,  60,   0, -60,  60,   0 /)	! "BO,3"

	basic_torsions(9,:)  = (/  30,  30, -60,  30,  30, -60 /)	! 1S5
	basic_torsions(10,:) = (/  60, -30, -30,  60, -30, -30 /)	! OS2
	basic_torsions(11,:) = (/  30, -60,  30,  30, -60,  30 /)	! 3S1
	basic_torsions(12,:) = (/ -30, -30,  60, -30, -30,  60 /)	! 5S1
	basic_torsions(13,:) = (/ -60,  30,  30, -60,  30,  30 /)	! 2SO
	basic_torsions(14,:) = (/ -30,  60, -30, -30,  60, -30 /)	! 1S3

	basic_torsions(15,:) = (/  45, -15,   0, -15,  45, -60 /)	! 1H2
	basic_torsions(16,:) = (/  60, -45,  15,   0,  15, -45 /)	! 3H2
	basic_torsions(17,:) = (/  45, -60,  45, -15,   0, -15 /)	! 3H4
	basic_torsions(18,:) = (/  15, -45,  60, -45,  15,   0 /)	! 5H4
	basic_torsions(19,:) = (/   0, -15,  45, -60,  45, -15 /)	! 5HO
	basic_torsions(20,:) = (/  15,   0,  15, -45,  60, -45 /)	! 1HO
	basic_torsions(21,:) = (/ -15,  45, -60,  45, -15,   0 /)	! 4H5
	basic_torsions(22,:) = (/   0,  15, -45,  60, -45,  15 /)	! OH5
	basic_torsions(23,:) = (/ -15,   0, -15,  45, -60,  45 /)	! OH1
	basic_torsions(24,:) = (/ -45,  15,   0,  15, -45,  60 /)	! 2H1
	basic_torsions(25,:) = (/ -60,  45, -15,   0, -15,  45 /)	! 2H3
	basic_torsions(26,:) = (/ -45,  60, -45,  15,   0,  15 /)	! 4H3

	basic_torsions(27,:) = (/  30,   0,   0, -30,  60, -60 /)	! 1E
	basic_torsions(28,:) = (/  60, -30,   0,   0,  30, -60 /)	! E2
	basic_torsions(29,:) = (/  60, -60,  30,   0,   0, -30 /)	! 3E
	basic_torsions(30,:) = (/  30, -60,  60, -30,   0,   0 /)	! E4
	basic_torsions(31,:) = (/   0, -30,  60, -60,  30,   0 /)	! 5E
	basic_torsions(32,:) = (/   0,   0,  30, -60,  60, -30 /)	! EO
	basic_torsions(33,:) = (/ -30,  60, -60,  30,   0,   0 /)	! 4E
	basic_torsions(34,:) = (/   0,  30, -60,  60, -30,   0 /)	! E5
	basic_torsions(35,:) = (/   0,   0, -30,  60, -60,  30 /)	! OE
	basic_torsions(36,:) = (/ -30,   0,   0,  30, -60,  60 /)	! E1
	basic_torsions(37,:) = (/ -60,  30,   0,   0, -30,  60 /)	! 2E
	basic_torsions(38,:) = (/ -60,  60, -30,   0,   0,  30 /)	! E3

	basicnames(:) = (/ "1C4","4C1", &
			&  "14B","B25","O3B","B14","25B","BO3", &
			&  "1S5","OS2","3S1","5S1","2SO","1S3", &
			&  "1H2","3H2","3H4","5H4","5HO","1HO","4H5","OH5","OH1","2H1","2H3","4H3", &
			&  "1E ","E2 ","3E ","E4 ","5E ","EO ","4E ","E5 ","OE ","E1 ","2E ","E3 " /)

	! Calculate norms
	do n=1,3
	  tempv = ideal_torsions(n,:) * ideal_torsions(n,:)
	  ideal_norms(n) = sqrt(sum(tempv))
	  tempv = redundant_torsions(n,:) * redundant_torsions(n,:)
	  redundant_norms(n) = sqrt(sum(tempv))
	end do
	do n=1,38
	  tempv = basic_torsions(n,:)*basic_torsions(n,:)
	  basic_norms(n) = sqrt(sum(tempv))
	end do

	end subroutine setupConformerData

	subroutine calculateConformer(torsions, d, phi, theta, verbose)
	implicit none
	real*8, intent(in) :: torsions(6)
	real*8, intent(out) :: d, phi, theta
	logical, intent(in) :: verbose
	real*8 :: lambda_ideal(3), lambda_redundant(3), lambda_basic(38), lambda_best(3), lambda_intermediate(3)
	real*8 :: remphi, tempf, r, bestr, bbtt
	real*8 :: tempv(6), normphi, lambda_c, lambda_b, lambda_t, norm_c, norm_b, norm_t
	integer :: best_basics(3), intermediate_basics(3), final_b, final_t, final_c, tempvi(6)
	integer :: n, m

	! 1) Work out projection coefficients of ideal (canonical) conformers and redundants by multiplying the
	!    torsion vector of the sugar by the relevant ideal/redundant torsion vector, and normalising
	do n=1,3
	  ! Ideals
	  tempv = torsions(:)*ideal_torsions(n,:)
	  lambda_ideal(n) = sum(tempv) / (ideal_norms(n)*ideal_norms(n))
	  ! Redundants
	  tempv = torsions(:)*redundant_torsions(n,:)
	  lambda_redundant(n) = sum(tempv) / (redundant_norms(n)*redundant_norms(n))
	end do
	if (verbose) write(6,"('Ideals (Redundants)',2x,3('[',f6.3,1x,a3,']',1x),3(f6.3,1x))") &
		& lambda_ideal(1),"1C4",lambda_ideal(2),"14B",lambda_ideal(3),"OS2",lambda_redundant(1:3)

	
	! 2) Determine best three conformers for description.
	! First, calculate projection coefficients of our sugar torsions on all 38 basic conformers
	do n=1,38
	  tempv = torsions(:)*basic_torsions(n,:)
	  lambda_basic(n) = sum(tempv) / (basic_norms(n)*basic_norms(n))
	end do
	! Now, determine best orthogonal boat-twistboat pair (boats: n=3,8, twistboats: 9-14)
	! 'Best' pair has highest 'R' coefficient.
	! Negative coefficients are ignored
	best_basics(BOAT) = 0
	best_basics(TWISTBOAT) = 0
	bestr = 0.0
	do n=14,3,-1
	  do m=14,3,-1
	    ! Determine if this pair are orthogonal
	    tempvi = basic_torsions(n,:)*basic_torsions(m,:)
	    if (sum(tempvi).eq.0) then
		! They are orthogonal as the projection coefficient is zero.
		! Now, work out normalized projection coefficient: (eq 7)
		! R(j) = lambda(i)norm(i) / lambda(j)norm(j)
	      r = ( lambda_basic(n)*basic_norms(n) ) / ( lambda_basic(m)*basic_norms(m) )
	!tempf = lambda_basic(BOAT) / lambda_basic(m)
	!write(6,"(a,a,1x,a,2x,4f8.3)") "Orthogonal pair: ",basicnames(n),basicnames(m),r,lambda_basic(n),lambda_basic(m),tempf
	      if ((lambda_basic(n).ge.0).and.(lambda_basic(m).ge.0).and.(r.ge.bestr)) then
		  ! The twistboat conformer may be in the 'BOAT' position, and vice versa...
		if (n.gt.8) then	! n is the twistboat, m is the boat
		  best_basics(BOAT) = m
		  best_basics(TWISTBOAT) = n
		else
		  best_basics(BOAT) = n
		  best_basics(TWISTBOAT) = m
		end if
	        bestr = r
	      end if
	    end if
	  end do
	end do
		

	! 3) Store projection coefficients of best boat/twistboat combination, and set the chair conformer
	lambda_best(BOAT) = lambda_basic(best_basics(BOAT))
	lambda_best(TWISTBOAT) = lambda_basic(best_basics(TWISTBOAT))
	! Select chair conformer with positive projection coefficient
	if (lambda_basic(1).ge.0) then
	  best_basics(CHAIR) = 1
	else
	  best_basics(CHAIR) = 2
	end if
	lambda_best(CHAIR) = lambda_basic(best_basics(CHAIR))
	if (verbose) write(6,"('	 Projected ',2x,3('[',f6.3,1x,a3,']',1x))") (lambda_best(n),basicnames(best_basics(n)),n=1,3)

	! 4) Check if we have an intermediate conformation
	! Copy the projection coefficients and conformers of the best basic trio
	lambda_intermediate = lambda_best
	intermediate_basics = best_basics
	tempf = lambda_intermediate(CHAIR) / lambda_intermediate(BOAT)
	if ((0.5.lt.tempf).and.(tempf.lt.2)) then
	  ! The chair component is now represented by 2*lambda(CHAIR) * 0.5(F(CHAIR) + F(BOAT))
	  ! So, we must find which basic conformation is equal to 0.5(F(CHAIR) + F(BOAT))...
	  tempvi = (basic_torsions(intermediate_basics(CHAIR),:) + basic_torsions(intermediate_basics(BOAT),:)) / 2
	  do n=1,38
	    r = sum(abs(tempvi(:) - basic_torsions(n,:)))
	!write(6,*) "intermediate search: ",basicnames(n),r
	    if (abs(r).lt.0.1) then
		! This is the conformation.
		intermediate_basics(CHAIR) = n
	!write(6,*) "Intermediate conformation found = ",basicnames(n),lambda_intermediate(CHAIR)
		exit
	    end if
	  end do
	  ! Set the new coefficients for boat and twistboat
	  lambda_intermediate(CHAIR) = 2 * lambda_best(CHAIR)
	  lambda_intermediate(BOAT) = lambda_best(CHAIR) - lambda_best(BOAT)
	  lambda_intermediate(TWISTBOAT) = lambda_best(TWISTBOAT)
	  ! Reorder boat and twistboat if necessary
	  if (lambda_intermediate(TWISTBOAT).gt.lambda_intermediate(BOAT)) then
	    r = lambda_intermediate(TWISTBOAT)
	    lambda_intermediate(TWISTBOAT) = lambda_intermediate(BOAT)
	    lambda_intermediate(BOAT) = r
	    n = intermediate_basics(TWISTBOAT)
	    intermediate_basics(TWISTBOAT) = intermediate_basics(BOAT)
	    intermediate_basics(BOAT) = n
	  end if
	  lambda_intermediate = abs(lambda_intermediate)
	  if (verbose) write(6,"('Needs intermediate ',2x,3('[',f6.3,1x,a3,']',1x))") (lambda_intermediate(n),basicnames(intermediate_basics(n)),n=1,3)
	end if
	  
	lambda_c = lambda_ideal(CHAIR)
	lambda_b = lambda_ideal(BOAT)
	lambda_t = lambda_ideal(TWISTBOAT)
	norm_c = ideal_norms(CHAIR)
	norm_b = ideal_norms(BOAT)
	norm_t = ideal_norms(TWISTBOAT)

	! 5) Convert these coefficients into the spherical polar coordinates
	! Notes:
	!  * Here the initial ratio is the calculation of phi TWISTBOAT / BOAT rather than the inverse
	!  * Explicit checks are made for 'perfect' chair conformers to prevent math errors
	!  * All calculations are made on the basic 1C4, 14B, and OS2 conformer projections

	! Meridian angle (phi) - Equation 12.
	! -- Check for exact chair conformers first
	if (dabs(lambda_b*norm_b).lt.1.0e-5) then
	  if ((lambda_t*norm_t).gt.0.0) then
	    phi = 90.0
	  else if ((lambda_t*norm_t).lt.0.0) then
	    phi = -90.0
	  else 
	    phi = 0.0
	  end if
	else 
	phi = ( lambda_t * norm_t ) / ( lambda_b * norm_b )
	phi = atan(phi)
	end if
	! Adjust phase of angle so that 14B is zero degrees.
	if (lambda_b.lt.0.0) then
	  phi = phi + pi
	else if ((lambda_t.lt.0.0).and.(lambda_b.gt.0.0)) then
	  phi = phi + 2.0*pi
	end if
	phi = phi * radcon
	
	! Calculate phi normalisation constant - Equation 14
	remphi = dabs(phi - (int(phi)/60)*60.0)
	normphi = sqrt( 2.0*(3600 + remphi*remphi + (60.0-remphi)*(60.0-remphi)))

	! Azimuthal angle (theta) - Equation 15.
	bbtt = (lambda_b*lambda_b) * (norm_b*norm_b) + (lambda_t*lambda_t) * (norm_t*norm_t)
	theta = (lambda_c*normphi) / sqrt(bbtt)
	theta = (0.5*pi - atan(theta))*radcon

	! Amplitude of sphere 'oblateness' - Equation 16
	d = sqrt((lambda_c * lambda_c) + (bbtt/(normphi*normphi)))
	if (verbose) then
	  write(6,"(a,f7.4,1x,a,f8.3,1x,a,f8.3)") "Amplitude =",d,"Phi =",phi,"Theta =",theta
	end if

	end subroutine calculateConformer

	end module conformerData


	program ringconf
	use dlprw; use utility; use parse; use conformerData; use PDensRW
	implicit none

	! Input / output files
	character*80 :: hisfile,dlpoutfile,torsfile,mapfile,altheaderfile,outputFile
	character*20 :: tempArg
	logical :: altheader = .false.
	integer :: iargc
	! Control
	integer :: n,m,m1,m2,p1,nframes,success,nargs
	integer :: framestodo = 0, nmapframes, startf = 0
	integer, allocatable :: molflags(:)
	logical :: molmap = .false.
	! Target Definition
	integer :: sp
	integer :: ringAtoms(6)
	real*8 :: testTorsions(6)
	! Calculation
	logical :: testCalculation = .false., torsionsFile = .false.
	real*8 :: torsions(6), d, phi, theta, x, y, z, mag, delta = 0.025
	integer :: i, j, k, l, moleculeCount = 0, grid, numAdded = 0
	type(PDens) :: conformers, conformersNorm
	logical :: addPoint

	nargs = iargc()
	if (nargs.lt.3) stop "Usage : ringconf <HISTORYfile> <OUTPUTfile> <sp> [-ring i j k l m n] [-test t1 t2 t3 t4 t5 t6] [-tors file]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	sp = getargi(3)
	
	n = 3
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,tempArg)
	  select case (tempArg)
            case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A,I4)") "Alternative header file supplied."
	      altheader = .TRUE.
	    case ("-molmap")
	      n = n + 1; call getarg(n,mapfile)
	      write(0,"(A,A)") "Using flags for species 1 molecules from : ",mapfile
	      open(unit=25,file=mapfile,form="formatted",status="old")
	      molmap = .true.
	    case ("-ring")
	      do m=1,6
		n = n + 1; ringAtoms(m) = getargi(n)
	      end do
	      write(0,"(a,6i4)") "Defined ring for target species : atoms ", ringAtoms
	    case ("-test")
	      do m=1,6
		n = n + 1; testTorsions(m) = getargr(n)
	      end do
	      write(0,"(a,6i4)") "Test torsions provided."
	      testCalculation = .true.
	    case ("-tors")
	      n = n + 1; call getarg(n,torsfile)
	      write(0,"(A,A)") "Reading torsions from : ", torsfile
	      open(unit=26,file=torsfile,form="formatted",status="old")
	      torsionsFile = .true.
	    case default
	      write(0,*) "Unrecognised command line option : ",tempArg
	      stop
	  end select
	end do

	! Perform initial setup
	call setupConformerData()

	grid = 1.5/delta
	conformers%axes = (/ delta,0.0d0,0.0d0,0.0d0,delta,0.0d0,0.0d0,0.0d0,delta /)
	conformers%origin = (/ -grid*delta,-grid*delta,-grid*delta /)
	if (.not.allocPDens(conformers,-grid,-grid,-grid,grid,grid,grid)) stop
	conformersNorm%axes = (/ delta,0.0d0,0.0d0,0.0d0,delta,0.0d0,0.0d0,0.0d0,delta /)
	conformersNorm%origin = (/ -grid*delta,-grid*delta,-grid*delta /)
	if (.not.allocPDens(conformersNorm,-grid,-grid,-grid,grid,grid,grid)) stop

	! If we are performing a test, set the torsions here, calculate the conformer, and exit
	if (testCalculation) then
	  write(6,"(a,6(f7.2,2x))") "Torsions: ", testTorsions
	  call calculateConformer(testTorsions, d, phi, theta, .true.)
	  call sphericalToCartesian(d, phi, theta, x, y, z)
	  write(6,"(a,3(f7.3,2x))") "Sphere mapping: ", x, y, z
	  write(6,*) ""

	  ! Add point to conformers grid
	  if (addPoint(conformers,x,y,z)) numAdded = numAdded + 1

	  ! Normalise to 1.0, and add to conformersNorm grid
	  mag = sqrt(x*x + y*y + z*z)
	  x = x / mag
	  y = y / mag
	  z = z / mag
	  if (addPoint(conformersNorm,x,y,z)) numAdded = numAdded + 1

	  goto 900

	end if

	! If we are reading from a torsions file, do it here and then exit
	if (torsionsFile) then
	  do while (readline(26))
	    if (nargsparsed.lt.6) then
	      write (6,"(a,i1,a)") "SKipping line with ", nargsparsed, " arguments on it"
	      continue
	    end if
	    do n=1,6
	      torsions(n) = argr(n)
	    end do
	    write(6,"(a,6(f7.2,2x))") "Torsions: ", torsions
	    call calculateConformer(torsions, d, phi, theta, .true.)
	    call sphericalToCartesian(d, phi, theta, x, y, z)
	    write(6,"(a,3(f7.3,2x))") "Sphere mapping: ", x, y, z
	    write(6,*) ""

	    ! Add point to conformers grid
	    if (addPoint(conformers,x,y,z)) numAdded = numAdded + 1

	    ! Normalise to 1.0, and add to conformersNorm grid
	    mag = sqrt(x*x + y*y + z*z)
	    x = x / mag
	    y = y / mag
	    z = z / mag
	    if (addPoint(conformersNorm,x,y,z)) numAdded = numAdded + 1
	  end do
	  goto 900
	end if

	! Open files
	if (outinfo(dlpoutfile,1).EQ.-1) goto 999

        ! Now, read in the history header so that we have cell()
	! If this fails then we may have a restarted trajectory. Continue, but only if
	! a header can be read in from the specified alternative history file..
	call openhis(hisfile,10)
        if (readheader().EQ.-1) then
	  if (altheader) then
	    write(0,*) "Restarted trajectory:"
	    close(dlpun_his)
	    call openhis(altheaderfile,10)
	    if (readheader().EQ.-1) goto 999
	    close(dlpun_his)
	    call openhis(hisfile,10)
	  else
	    goto 999
	  end if
	end if

	if (molmap) allocate(molflags(s_nmols(sp)))

	! XXXX
	! XXXX Main routine....
	! XXXX
	nframes=0
	nmapframes=0
	! If we're using a mapfile, read it here
101     if (molmap) then
	  read(25,"(I5,1x,40(I2,1x))",end=118,err=117) m1,(molflags(m2),m2=1,s_nmols(sp))
	  ! if (m1.NE.nframes) stop "Mapfile missed a frame number!"
	end if

102	success=readframe()
	if (success.eq.1) goto 800  ! End of file encountered....
	if (success.lt.0) goto 800  ! File error....

	nframes=nframes+1
	if (mod(nframes,100).EQ.0) then
	  if (molmap) then
	    write(0,"(I6,2x,'(',I6,',',I6,')')") nframes,nmapframes
	  else
	    write(0,"(i6)") nframes
	  end if
	end if

	if (nframes.lt.startf) goto 102

	if ((molmap).and.(m1.ne.nframes)) goto 102
	if (molmap) nmapframes = nmapframes + 1

	! Loop over molecules of the target type
	p1 = s_start(sp) - s_natoms(sp) - 1
	do m1 = 1,s_nmols(sp)
	  p1 = p1 + s_natoms(sp)

	  ! MolMap check
	  if (molmap.and.(molflags(m1).eq.0)) cycle

	  moleculeCount = moleculeCount + 1

	  ! Calculate torsions for this molecule
	  do n=1,6
	    i = ringAtoms(mod(n-1,6)+1)
	    j = ringAtoms(mod(n,6)+1)
	    k = ringAtoms(mod(n+1,6)+1)
	    l = ringAtoms(mod(n+2,6)+1)
	    torsions(n) = calculateTorsion(i, j, k, l, p1)
	  end do

	  ! Get spherical and cartesian coordinates
	  call calculateConformer(torsions, d, phi, theta, .false.) 
	  call sphericalToCartesian(d, phi, theta, x, y, z)

	  ! Add point to conformers grid
	  if (addPoint(conformers,x,y,z)) numAdded = numAdded + 1

	  ! Normalise to 1.0, and add to conformersNorm grid
	  mag = sqrt(x*x + y*y + z*z)
	  x = x / mag
	  y = y / mag
	  z = z / mag
	  if (addPoint(conformersNorm,x,y,z)) numAdded = numAdded + 1

	end do

	if (nframes.EQ.framestodo) goto 800

	! Next frame
	goto 101

117     write(0,*) "Error reading mapfile."
	write(0,*) "Selected ",nmapframes," from trajectory before error."
	goto 801
118     write(0,*) "Reached end of mapfile."
	write(0,*) "Selected ",nmapframes," from trajectory."
	goto 801
799	write(0,*) "Error in data files."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",framestodo,")"
	goto 801
800	write(0,*) "Framestodo was fulfilled."
801	write(0,*) ""

900	write(0,*) "Writing output files..."

	conformers%grid = conformers%grid / (numAdded)
	outputFile = outputFileName(hisfile, "test", "pdens")
	if (.not.savePDens(outputFile, conformers)) write(0,*) "Error saving pdens file."

	conformersNorm%grid = conformersNorm%grid / (numAdded)
	outputFile = outputFileName(hisfile, "test", "norm.pdens")
	if (.not.savePDens(outputFile, conformersNorm)) write(0,*) "Error saving pdens file."

999	write(0,*) "Finished."
	end program ringconf

	subroutine sphericalToCartesian(d,phi,theta,x,y,z)
	use conformerData
	implicit none
	real*8, intent(in) :: d, phi, theta
	real*8, intent(out) :: x, y, z
	real*8 :: thetaRad, phiRad

	! Theta = polar angle / inclination (angle deviation from O-Z)
	! Phi = azimuthal angle (angle rotation in XY plane)
	! Assumes that phi and theta are in degrees
	! Note: Signs of x and y coordinates are switched to give B1,4 along positive X, and 2SO along positive Y
	thetaRad = theta / radcon
	phiRad = phi / radcon
	x = -d * sin(thetaRad) * cos(phiRad);
	y = -d * sin(thetaRad) * sin(phiRad);
	z = d * cos(thetaRad);

	end subroutine sphericalToCartesian

	logical function addPoint(grid,x,y,z)
	use PDensRW
	implicit none
	type(PDens), intent(inout) :: grid
	real*8, intent(in) :: x, y, z
	integer :: nx, ny, nz

	! Work out grid coordinates, assuming that we have orthogonal voxel axes
	nx = x / grid%axes(1)
	ny = y / grid%axes(5)
	nz = z / grid%axes(9)
	addPoint = .false.
	if ((nx.lt.grid%gridMin(1)).or.(nx.gt.grid%gridMax(1))) return
	if ((ny.lt.grid%gridMin(2)).or.(ny.gt.grid%gridMax(2))) return
	if ((nz.lt.grid%gridMin(3)).or.(nz.gt.grid%gridMax(3))) return

	grid%grid(nx,ny,nz) = grid%grid(nx,ny,nz) + 1
	addPoint = .true.
	end function addPoint
