	! ################################################################
	! geom - analyses the average geometry in a history or config file
	! ################################################################

	program geom
	use utility; use dlprw; use parse
	implicit none
	integer, parameter :: MAXDATA = 100
	real*8, dimension(:,:), allocatable :: gadata, gbdata, gtdata
	real*8, allocatable :: calctors(:)
	real*8  :: tx,ty,tz,vecji(3),vecjk(3),veckj(3),veckl(3),xp1(3),xp2(3),mag1,mag2,ktx,kty,ktz
        real*8  :: dp,angle,radcon,minimum,maximum,avg,sd,total, rij
	character*80 :: outfile,hisfile,temp,altheaderfile
	integer :: nargs,n,m,i,j,k,l,status, nframes, framestodo, sp, a1, ngeom, num
	integer :: nbonds, nangles, ntorsions, success, framestodiscard = 0
	integer :: bonds(MAXDATA,2),angles(MAXDATA,3),torsions(MAXDATA,4)
	integer :: iargc
	logical :: configfile = .false., altheader = .false.
	real*8 :: warndist = 100.0, warnangle = 200.0
1       FORMAT (A8/,3f20.14)
2       FORMAT (3F20.14)
3	FORMAT (2I10)

	radcon = 57.29577951d0

	nargs = iargc()
	if (nargs.lt.5) stop "Usage : geom <HISfile|CONFIG> <OUTfile> <sp> <nframes> [-bond i j] [-angle i j k] [-torsion i j k l] [-data file] [-discard n] [-warndist maxdist] [-warnangle maxangle]"
	call getarg(1,hisfile)
	call getarg(2,outfile)
	sp = getargi(3)
	framestodo = getargi(4)

	nbonds = 0
	nangles = 0
	ntorsions = 0

	n = 4
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
	    case ("-discard")
              n = n + 1; framestodiscard = getargi(n)
            case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
	      altheader = .TRUE.
	    case ("-warndist")
              n = n + 1; warndist = getargr(n)
	    case ("-warnangle")
              n = n + 1; warnangle = getargr(n)
            case ("-bond")
	      nbonds = nbonds + 1
	      if (nbonds.gt.MAXDATA) stop "Bond array MAXDATA exceeded. Increase and recompile."
              n = n + 1; bonds(nbonds,1) = getargi(n)
              n = n + 1; bonds(nbonds,2) = getargi(n)
            case ("-angle")
	      nangles = nangles + 1
	      if (nangles.gt.MAXDATA) stop "Angle array MAXDATA exceeded. Increase and recompile."
              n = n + 1; angles(nangles,1) = getargi(n)
              n = n + 1; angles(nangles,2) = getargi(n)
              n = n + 1; angles(nangles,3) = getargi(n)
            case ("-torsion")
	      ntorsions = ntorsions + 1
	      if (ntorsions.gt.MAXDATA) stop "Torsion array MAXDATA exceeded. Increase and recompile."
              n = n + 1; torsions(ntorsions,1) = getargi(n)
              n = n + 1; torsions(ntorsions,2) = getargi(n)
              n = n + 1; torsions(ntorsions,3) = getargi(n)
              n = n + 1; torsions(ntorsions,4) = getargi(n)
	    case ("-data")
	      n = n + 1; call getarg(n,temp)
	      open(unit=20, file=temp, form='formatted', status='old', err=999)
	      do while (readline(20))
		select case (arg(1))
		  case ("bond")
		    nbonds = nbonds + 1
		    if (nbonds.gt.MAXDATA) stop "Bond array MAXDATA exceeded. Increase and recompile."
		    bonds(nbonds,1) = argi(2); bonds(nbonds,2) = argi(3)
		  case ("angle")
		    nangles = nangles + 1
		    if (nangles.gt.MAXDATA) stop "Angle array MAXDATA exceeded. Increase and recompile."
		    angles(nangles,1) = argi(2); angles(nangles,2) = argi(3); angles(nangles,3) = argi(4)
		  case ("torsion")
		    ntorsions = ntorsions + 1
		    if (ntorsions.gt.MAXDATA) stop "Torsion array MAXDATA exceeded. Increase and recompile."
		    torsions(ntorsions,1) = argi(2); torsions(ntorsions,2) = argi(3); torsions(ntorsions,3) = argi(4); torsions(ntorsions,4) = argi(5)
		  case default
		    write(0,"(a,a)") "Unrecognised geometry type:",arg(1)
		    stop
		end select
	      end do
	      close(20)
	    case default
	      write(0,"(a,a)") "Unrecognised command line option:",temp
	      stop
	  end select
	end do

	! Open and check the files...
	if (outinfo(outfile,1).EQ.-1) stop "Bad OUTPUT file."


	if (nbonds.gt.0) then
	  allocate(gbdata(nbonds,framestodo*s_nmols(sp)))
	  write(6,*) nbonds,"bonds defined."
	  gbdata = 0.0
	end if
	if (nangles.gt.0) then
	  allocate(gadata(nangles,framestodo*s_nmols(sp)))
	  write(6,*) nangles,"angles defined."
	  gadata = 0.0
	end if
	if (ntorsions.gt.0) then
	  allocate(gtdata(ntorsions,framestodo*s_nmols(sp)))
	  write(6,*) ntorsions,"torsions defined."
	  gtdata = 0.0
	end if

	! Open the history file or config file
	if ((index(hisfile,"CONFIG").ne.0).or.(index(hisfile,"REVCON").ne.0)) then
	  configfile = .true.
	  framestodiscard = 0
	  framestodo = 1
	  call readconfig(hisfile)
	else
	  call openhis(hisfile,15)
          if (readheader().EQ.-1) then
	    if (altheader) then
	      write(0,*) "Restarted trajectory:"
	      close(dlpun_his)
	      call openhis(altheaderfile,15)
	      if (readheader().EQ.-1) stop "Error reading header / history file."
	      close(dlpun_his)
	      call openhis(hisfile,15)
	    else
	      stop "Error reading header / history file."
	    end if
	  end if
	end if

	nframes = 0
	num = 0

10	if (configfile) then
	  call readconfig(hisfile)
	  do n=1,natms
	    xpos(n) = cfgxyz(n,1)
	    ypos(n) = cfgxyz(n,2)
	    zpos(n) = cfgxyz(n,3)
	  end do
	  cell = cfgcell
	  success = 0
	else
	  success = readframe()
	end if

	if (success.eq.1) then
	  write(0,*) "End of HISTORY file found."
	  goto 19
	end if
	if (success.ne.0) then
	  write(0,*) "Error reading history file. Error code:", success
	  goto 19
	end if
	if (framestodiscard.gt.0) then
	  framestodiscard = framestodiscard -1
	  goto 10
	end if

	nframes = nframes + 1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	a1 = s_start(sp) - 1
	do m=1,s_nmols(sp)
	  num = num + 1

	  ! Calculate bonds
	  do n=1,nbonds
	    i = bonds(n,1) + a1
	    j = bonds(n,2) + a1
	    call pbc(xpos(i),ypos(i),zpos(i),xpos(j),ypos(j),zpos(j),tx,ty,tz)
	    rij = sqrt( (tx-xpos(j))*(tx-xpos(j)) + (ty-ypos(j))*(ty-ypos(j)) + (tz-zpos(j))*(tz-zpos(j)) )
	    gbdata(n,num) = rij
	    if (rij.gt.warndist) then
	      write(6,"(a,i4,'/',i6,'/',i6,'/',f12.6)") "Big bond found : mol/i/j/rij = ", m, i, j, rij
	      write(6,"(9f12.6)") xpos(i),ypos(i),zpos(i),tx,ty,tz,xpos(j),ypos(j),zpos(j)
	    end if
	  end do

	  ! Calculate angles....
	  do n=1,nangles
	    angle = calculateAngle(angles(n,1),angles(n,2),angles(n,3),a1)
	    if (angle.gt.warnangle) write(6,"(a,i4,'/',i6,'/',i6,'/',i6,'/',f12.6)") "Big angle found : mol/i/j/k/rij = ", m, i, j, k, angle
	    gadata(n,num) = angle
	  end do

	  ! Calculate torsions....
	  do n=1,ntorsions
	    gtdata(n,num) = calculateTorsion(torsions(n,1),torsions(n,2),torsions(n,3),torsions(n,4),a1)
          end do

	  a1 = a1 + s_natoms(sp)
	end do

	if (nframes.EQ.framestodo) goto 19
	goto 10


19	write(0,"(A)") "Analysing data...."
	close(14)
20	FORMAT ('Bond    ',I2,' : Atoms ',6x,2I3,4(2x,f12.6))
21	FORMAT ('Angle   ',I2,' : Atoms ',3x,3I3,4(2x,f12.6))
22	FORMAT ('Torsion ',I2,' : Atoms ',4I3,4(2x,f12.6))
23	FORMAT ('                                         Min           Max           Avg          S.D.')

	! Calculate statistics
	write(6,23)
	do n=1,nbonds
	  minimum = minval(gbdata(n,1:num))
	  maximum = maxval(gbdata(n,1:num))
	  avg = sum(gbdata(n,1:num)) / real(num)
	  ! Now for S.D.
	  total = 0.0d0
	  do m=1,num
	    total = total + (gbdata(n,m) - avg)**2
	  end do
	  sd = SQRT( total / real(num) )
	  write(6,20) n,(bonds(n,m),m=1,2),minimum,maximum,avg,sd
	end do

	do n=1,nangles
	  minimum = minval(gadata(n,1:num))
	  maximum = maxval(gadata(n,1:num))
	  avg = sum(gadata(n,1:num)) / real(num)
          ! Now for S.D.
          total = 0.0d0
          do m=1,num
            total = total + (gadata(n,m) - avg)**2
          end do
          sd = SQRT( total / real(num) )
          write(6,21) n,(angles(n,m),m=1,3),minimum,maximum,avg,sd
        end do

	do n=1,ntorsions
	  minimum = minval(gtdata(n,1:num))
	  maximum = maxval(gtdata(n,1:num))
	  avg = sum(gtdata(n,1:num)) / real(num)
          ! Now for S.D.
          total = 0.0d0
          do m=1,num
            total = total + (gtdata(n,m) - avg)**2
          end do
          sd = SQRT( total / real(num) )
          write(6,22) n,(torsions(n,m),m=1,4),minimum,maximum,avg,sd
        end do

999	stop
	end program geom
	
