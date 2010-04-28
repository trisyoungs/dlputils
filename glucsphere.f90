	! glucsphere
	! Calculates the spherical polar coordinate description of the sugar's conformation, using the method
	! of Berces et al.

	program glucsphere
	implicit none

	integer, parameter :: CHAIR=1, BOAT=2, TWISTBOAT=3
	real*8, parameter :: radcon = 57.29577951d0
	real*8, parameter :: pi = 3.14159265358979d0
	character*80 :: torsfile,flagfile
	character*20 :: temparg
	character*3 :: basicnames(38)
	integer :: n,nframes,success,nargs,m,cl,gluc,m1,m2
	integer :: framestodo = 0, ntorsions, nglucose, nmapframes, startf = 0, tempvi(6)
	integer :: best_basics(3), intermediate_basics(3), final_b, final_t, final_c
	integer, allocatable :: molflags(:)
	integer :: iargc
	logical :: molmap
	real*8, allocatable :: torsions(:,:)
	integer :: ideal_torsions(3,6), redundant_torsions(3,6), basic_torsions(38,6)
	real*8 :: ideal_norms(3), redundant_norms(3), basic_norms(38)
	real*8 :: lambda_ideal(3), lambda_redundant(3), lambda_basic(38), lambda_best(3), lambda_intermediate(3)
	real*8 :: tempv(6), normphi, lambda_c, lambda_b, lambda_t
	real*8 :: phi, theta, d, remphi, tempf, r, bestr, offsets(38)

	nargs = iargc()
	if (nargs.ne.2) stop "Usage : glucsphere <torsions.dat> <ngluc>"
	call getarg(1,torsfile)
	call getarg(2,temparg)
	read(temparg,"(I5)") nglucose
	
	n = 4
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temparg)
	  select case (temparg)
	    case ("-molmap")
	      n = n + 1; call getarg(n,flagfile)
	      write(0,"(A,A)") "Using flags for species 1 molecules from : ",flagfile
	      open(unit=25,file=flagfile,form="formatted",status="old")
	      molmap = .TRUE.
	    case default
	      write(0,*) "Unrecognised command line option : ",temparg
	      stop
	  end select
	end do

	! Open files
	open(unit=19,file=torsfile,form="formatted",status="old")
	open(unit=20,file="sphericalcoords.dat",form="formatted",status="replace")

	! Set some numbers
	ntorsions = nglucose * 6

	! Set ideal, redundancy, and basic torsion vectors.
	! Redundant torsions R4, R5, and R6 map to redundant_torsions(1-3,:)
	ideal_torsions(CHAIR,:)     = (/ 60, -60, 60, -60, 60, -60 /)	! 1C4
	ideal_torsions(BOAT,:)      = (/ 0, 60, -60, 0, 60, -60 /)	! 14B
	ideal_torsions(TWISTBOAT,:) = (/ 60, -30, -30, 60, -30, -30 /)	! OS2
	redundant_torsions(1,:)     = (/ 60, 30, -30, -60, -30, 30 /)
	redundant_torsions(2,:)     = (/ 0, 60, 60, 0, -60, -60 /)
	redundant_torsions(3,:)     = (/ 60, 60, 60, 60, 60, 60 /)

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

	offsets(:) = (/ 90.0, 180.0,    0.0, 60.0, 120.0, 180.0, 240.0, 300.0,    30.0, 90.0, 150.0, 210.0, 270.0, 330.0, &
			& 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

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

	if (molmap) allocate(molflags(nglucose))
	allocate(torsions(nglucose,6))

	! XXXX
	! XXXX Main routine....
	! XXXX
	nframes=0
	nmapframes=0
	! If we're using a mapfile, read it here
101     if (molmap) then
	  read(25,"(I5,1x,40(I2,1x))",end=118,err=117) m1,(molflags(m2),m2=1,nglucose)
	  ! if (m1.NE.nframes) stop "Mapfile missed a frame number!"
	end if
102	read(19,"(11x,100(f5.1,1x))",end=801,err=799) ((torsions(n,m),m=1,6),n=1,nglucose)

	nframes=nframes+1
	if (mod(nframes,100).EQ.0) then
	  if (molmap) then
	    write(0,"(I6,2x,'(',I6,',',I6,')')") nframes,nmapframes
	  else
	    write(0,"(i6)") nframes
	  end if
	end if

	if (nframes.LT.startf) goto 102

	if ((molmap).and.(m1.ne.nframes)) goto 102
	if (molmap) nmapframes = nmapframes + 1

	! CALCULATION STARTS HERE

	select case (nframes)
	  case (1)
		torsions(1,:) = (/ 42.7672, -42.7123, 54.1549, -66.6666, 65.8853, -52.921 /)	! Whitfield [1] (Ref 30)
	  case (2)
		torsions(1,:) = (/ 61.2798, -24.1404, -36.4885, 70.8247, -34.0682, -31.8185 /)	! Chamberlain [2] (Ref 31)
	  case (3)
		torsions(1,:) = (/ -73.3374, 38.4315, 0.8983, -0.5636, -39.0175, 75.5659 /)	! Urbanczyk [3] (Ref 32)
	  case (4)
		torsions(1,:) = (/ 45.0, -15.0, 0.0, -15.0, 45.0, -60.0 /)			! 1H2
	  case (5)
		torsions(1,:) = (/60.0,-60.0,   0.0,60.0,-60.0,   0.0 /)
	  case (6)
		torsions(1,:)  = (/   0, -60,  60,   0, -60,  60 /)	! "B1,4"
	  case (7)
		torsions(1,:) = (/ -60,  60, -60,  60, -60,  60 /)	! 4C1
	  case (8)
		torsions(1,:) = (/ -30, -30,  60, -30, -30,  60 /)	! 5S1
	end select

	do gluc=1,nglucose

	  write(0,"(A,6f6.2)") "Torsions = ",torsions(gluc,:)
	  ! 1) Work out projection coefficients of ideal (canonical) conformers and redundants by multiplying the
	  !    torsion vector of the sugar by the relevant ideal/redundant torsion vector, and normalising
	  do n=1,3
	    ! Ideals
	    tempv = torsions(gluc,:)*ideal_torsions(n,:)
	    lambda_ideal(n) = sum(tempv) / (ideal_norms(n)*ideal_norms(n))
	    ! Redundants
	    tempv = torsions(gluc,:)*redundant_torsions(n,:)
	    lambda_redundant(n) = sum(tempv) / (redundant_norms(n)*redundant_norms(n))
	  end do
	  write(0,"('Ideals (Redundants)',2x,3('[',f6.3,1x,a3,']',1x),3(f6.3,1x))") &
		& lambda_ideal(1),"1C4",lambda_ideal(2),"14B",lambda_ideal(3),"OS2",lambda_redundant(1:3)

	
	  ! 2) Determine best three conformers for description.
	  ! First, calculate projection coefficients of our sugar torsions on all 38 basic conformers
	  do n=1,38
	    tempv = torsions(gluc,:)*basic_torsions(n,:)
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
	!write(0,"(a,a,1x,a,2x,4f8.3)") "Orthogonal pair: ",basicnames(n),basicnames(m),r,lambda_basic(n),lambda_basic(m),tempf
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
	  write(0,"('	 Projected ',2x,3('[',f6.3,1x,a3,']',1x))") (lambda_best(n),basicnames(best_basics(n)),n=1,3)

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
	!write(0,*) "intermediate search: ",basicnames(n),r
	      if (abs(r).lt.0.1) then
		! This is the conformation.
		intermediate_basics(CHAIR) = n
	!write(0,*) "Intermediate conformation found = ",basicnames(n),lambda_intermediate(CHAIR)
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
	    write(0,"('Needs intermediate ',2x,3('[',f6.3,1x,a3,']',1x))") (lambda_intermediate(n),basicnames(intermediate_basics(n)),n=1,3)
	    ! Store intermediate data in best_basics and lambda_best so it is used subsequently...
	!    best_basics = intermediate_basics
	!    lambda_best = lambda_intermediate
	  end if
	    
	  final_c = best_basics(CHAIR)
	  final_b = best_basics(BOAT)
	  final_t = best_basics(TWISTBOAT)
	  lambda_c = lambda_best(CHAIR)
	  lambda_b = lambda_best(BOAT)
	  lambda_t = lambda_best(TWISTBOAT)

	  ! 5) Convert these coefficients into the spherical polar coordinats
	  ! Meridian angle (phi) - Equation 12.
	  phi = ( lambda_b * lambda_b ) / ( lambda_t * lambda_t )
	  phi = atan(phi) * radcon
	  ! Adjust angle from the 'local' frame of the best boat-twistboat pair into the global frame where 14B is 0deg.
	  phi = phi + offsets(final_t)
	  ! Calculate phi normalisation constant
	  remphi = phi - (int(phi)/60)*60.0
	  normphi = sqrt( 2.0*3600 + 2.0*(remphi*remphi) + 2*((60.0-remphi)*(60.0-remphi)))
	  ! Azimuthal angle (theta) - Equation 15.
	  tempf = lambda_b*lambda_b * basic_norms(final_b)*basic_norms(final_b)
	  tempf = tempf + lambda_t*lambda_t * basic_norms(final_t)*basic_norms(final_t)
	  theta = ( lambda_c*normphi ) / sqrt(tempf)
	  theta = offsets(final_c) - atan(theta)*radcon
	  ! Amplitude of sphere 'oblateness'
	!  tempf = lambda_intermediate(BOAT)*lambda_intermediate(BOAT) * basic_norms(intermediate_basics(BOAT))*basic_norms(intermediate_basics(BOAT))
	  !tempf = tempf + lambda_intermediate(TWISTBOAT)*lambda_intermediate(TWISTBOAT) * basic_norms(intermediate_basics(TWISTBOAT))*basic_norms(intermediate_basics(TWISTBOAT))
	  !d = sqrt( lambda_intermediate(CHAIR)*lambda_intermediate(CHAIR) + tempf/(normphi*normphi) )
	d = 0.0
	  write(0,"(a,f7.4,1x,a,f8.3,1x,a,f8.3)") "Amplitude =",d,"Phi =",phi,"Theta =",theta
	  write(0,*) ""
	
	end do	! End loop over glucose molecules

	if (nframes.eq.3) stop
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

	close(20)
	close(19)

999	write(0,*) "Finished."
	end program glucsphere
