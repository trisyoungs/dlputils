!	** point.f90 **
!	Program to calculate the vector of molecules (pointing vector) in a system

	program point
	use dlprw; use utility
	implicit none
	!integer, parameter :: ngrid=100
	integer :: ngrid=100
	integer :: maxframes,vec(3),baselen,nargs,success,n,i,targetsp,startf,endf
	integer :: nframes, nframesused, xbin, ybin, zbin, zspecies
	integer :: anglebin1, anglebin2, framestodiscard=-1, nskipped
	character*80 :: hisfile,outfile,basename,resfile,temp,flagfile,altheaderfile
	character*8 :: discard
	logical :: altheader = .FALSE., usecom = .FALSE., doangles =.FALSE.
	real*8 :: histo(-100:100) !acc(-ngrid:ngrid,-ngrid:ngrid), anglehist(-ngrid:ngrid,-ngrid:ngrid,2,0:180), map(-ngrid:ngrid,-ngrid:ngrid,5)
	real*8,allocatable :: acc(:,:), anglehist(:,:,:,:), map(:,:,:), legmap(:,:)
	real*8 :: vec1(3), vec2(3), pvec(3), mag, tx, ty, tz, pvecmult, theta
	real*8 :: minz, maxz, xdelta, ydelta, patchz, globalacc, globalsum, globalsumpos, zorigin, vel
	real*8, parameter :: radcon = 57.29577951d0
	integer :: iargc

	write(0,*) "*** pdens"

	nargs = iargc()
	if (nargs.LT.2) then
	  write(*,"(a)") "Usage: point <HIStory file> <OUTput file> [-options]"
	  write(*,"(a)") "        [-target sp]            Set species sp to be the target"
	  write(*,"(a)") "        [-vector x1 x2 y1]      Atoms to use for vector calculation in target species"
	  write(*,"(a)") "        [-zspecies sp]          Use specified species average Z as first distance point"
	  write(*,"(a)") "        [-zorigin z]            Use specified Z as first distance point"
	  write(*,"(a)") "        [-minz r]               Minimum deviation from first distance point"
	  write(*,"(a)") "        [-maxz r]               Maximum deviation from first distance point"
	  write(*,"(a)") "        [-header file]          Use specified file to get header"
	  write(*,"(a)") "        [-grid n]               Use specified number of gridpoints in X/Y"
	  write(*,"(a)") "        [-reversevec]           Reverse sign of vector"
	  write(*,"(a)") "        [-discard n]               Skip n frames at start of trajectory"
	  write(*,"(a)") "        [-com]                  Use centre of mass for z-coordinate check rather than vector centre"
	  write(*,"(a)") "        [-angles]               Calculate surface-molecule distribution map"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
	end if

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  IF (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 50
	  endIF
	end do
50	if (baselen.eq.-1) then
	  basename="pointresults."
	  baselen=14
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=20,file=basename(1:baselen)//"out",status='replace',form='formatted')

	write(0,"(A,A)") "History file : ",hisfile
	write(20,"(A,A)") "History file : ",hisfile
	!call openhis(hisfile,n)
	write(0,"(A,A)") " Output file : ",outfile
	write(20,"(A,A)") " Output file : ",outfile
	if (outinfo(outfile,1).eq.-1) goto 798

	! Set some variable defaults before we read in any command-line arguments
	pvecmult = -1.0
	vec(1) = 1
	vec(2) = 2
	vec(3) = 3
	minz = 0.0
	maxz = 100.0
	zorigin = 0.0
	zspecies = 0
	targetsp = 1
	histo = 0.0
	globalacc = 0.0
	globalsum = 0.0

	n = 2
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-vector") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") vec(1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") vec(2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") vec(3)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Vector for species ",targetsp," calculated from: A=",vec(1),"->", &
	        & vec(2),", B=",vec(3)
	      write(20,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Vector for species ",targetsp," calculated from: A=",vec(1),"->", &
	        & vec(2),", B=",vec(3)
	    case ("-target")
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") targetsp
	      write(0,"(A,I4)") "Target species is = ",targetsp
	      write(20,"(A,I4)") "Target species is = ",targetsp
	    case ("-maxz")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") maxz
	      write(0,"(A,f8.4)") "Maximum z deviation = ",maxz
	      write(20,"(A,f8.4)") "Maximum z deviation = ",maxz
	    case ("-minz")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") minz
	      write(0,"(A,f8.4)") "Minimum z deviation = ",minz
	      write(20,"(A,f8.4)") "Minimum z deviation = ",minz
	    case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
              write(20,"(A)") "Alternative header file supplied."
              altheader = .TRUE.
	    case ("-zspecies")
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") zspecies
	      write(0,"(A,I4)") "Z-basis species is = ",zspecies
	      write(20,"(A,I4)") "Z-basis species is = ",zspecies
	    case ("-zorigin")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") zorigin
	      write(0,"(A,f8.4)") "Z-basis origin is = ",zorigin
	      write(20,"(A,f8.4)") "Z-basis origin is = ",zorigin
	    case ("-grid")
	      n = n + 1; call getarg(n,temp); read(temp,"(i5)") ngrid
	      write(0,"(A,i6)") "Number of grid points in each X/Y direction = ",ngrid
	      write(20,"(A,i6)") "Number of grid points in each X/Y direction = ",ngrid
	    case ("-reversevec")
	      write(0,"(A)") "Sign of pointing vector will be negated."
	      write(20,"(A)") "Sign of pointing vector will be negated."
	      pvecmult = 1.0
	    case ("-com")
	      write(0,"(A)") "Molecule COM will be used for z-coordinate."
	      write(20,"(A)") "Molecule COM will be used for z-coordinate."
	      usecom = .TRUE.
	    case ("-angles")
	      write(0,"(A)") "Surface-molecule angle distributions will be calculated."
	      write(20,"(A)") "Surface-molecule angle distributions will be calculated."
	      doangles = .TRUE.
	    case ("-discard")
	      n = n + 1; call getarg(n,temp); read(temp,"(i5)") framestodiscard
	      write(0,"(A,i6)") "Frames to discard at start of file:",framestodiscard
	      write(20,"(A,i6)") "Frames to discard at start of file:",framestodiscard
	    case default
	      write(0,*) "Unrecognised command line option : ",temp
	      stop
	  end select
	end do
	
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

	! Allocate arrays
	allocate(acc(-ngrid:ngrid,-ngrid:ngrid))
	if (doangles) allocate(anglehist(-ngrid:ngrid,-ngrid:ngrid,2,0:180))
	allocate(map(-ngrid:ngrid,-ngrid:ngrid,5))
	if (doangles) anglehist = 0.0
	map = 0.0
	acc = 0.0

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes = 0
	nskipped = 0
	nframesused = 0
	xdelta = cell(1) / real(2*ngrid+1)
	ydelta = cell(5) / real(2*ngrid+1)

102	success=readframe()
	IF (success.eq.1) goto 801  ! End of file encountered....
	IF (success.eq.-1) goto 119  ! File error....
	nskipped = nskipped + 1
	if (nskipped.le.framestodiscard) goto 102
	nframes=nframes+1
	if (MOD(nframes,100).eq.0) write(0,*) nframes
	if (nframes.LT.startf) goto 102

	nframesused = nframesused + 1

	! Calculate COM if needed
	if (usecom) call calc_com()

	! Create patch position
	if ((zspecies.ne.0).and.(nframesused.eq.1)) then
	  write(20,*) "Writing molecule info..."
	  write(20,*) s_nmols(zspecies)
	  write(20,*) "Patch"
	  i = s_start(zspecies)
	  do n=1,s_nmols(zspecies)
	    ! Get first atom of point A
	    tx = xpos(i+n-1)
	    ty = ypos(i+n-1)
	    tz = zpos(i+n-1)
	    ! Fold into -0.5->0.5 (mostly for first frame in DLP hisfile)
	    if (tx.lt.(-0.5*cell(1))) tx = tx + cell(1)
	    if (tx.gt.(0.5*cell(1))) tx = tx - cell(1)
	    if (ty.lt.(-0.5*cell(5))) ty = ty + cell(5)
	    if (ty.gt.(0.5*cell(5))) ty = ty - cell(5)
	    if (tz.lt.(-0.5*cell(9))) tz = tz + cell(9)
	    if (tz.gt.(0.5*cell(9))) tz = tz - cell(9)
	    write(20,*) "C",tx,ty,0.0
	  end do
	endif
	if (zspecies.eq.0) then
	  patchz = zorigin
	else
	  patchz = zpos(s_start(zspecies))
	end if
	if (patchz.lt.(-0.5*cell(9))) patchz = patchz + cell(9)
	if (patchz.gt.(0.5*cell(9))) patchz = patchz - cell(9)

	i = s_start(targetsp)
	do n=1,s_nmols(targetsp)

	  ! Get mim unit vector position1 from first atom with third atom (in vec1)
	  call pbc(xpos(i+vec(1)-1),ypos(i+vec(1)-1),zpos(i+vec(1)-1), &
		& xpos(i+vec(3)-1),ypos(i+vec(3)-1),zpos(i+vec(3)-1), vec1(1), vec1(2), vec1(3))
	  vec1(1) = vec1(1) - xpos(i+vec(3)-1)
	  vec1(2) = vec1(2) - ypos(i+vec(3)-1)
	  vec1(3) = vec1(3) - zpos(i+vec(3)-1)
	  mag = sqrt(vec1(1)*vec1(1) + vec1(2)*vec1(2) + vec1(3)*vec1(3))
	  vec1(1) = vec1(1) / mag
	  vec1(2) = vec1(2) / mag
	  vec1(3) = vec1(3) / mag
	  
	  ! Get MIM unit vector position from second atom with third atom (in vec2)
	  call pbc(xpos(i+vec(2)-1),ypos(i+vec(2)-1),zpos(i+vec(2)-1), &
		& xpos(i+vec(3)-1),ypos(i+vec(3)-1),zpos(i+vec(3)-1), vec2(1), vec2(2), vec2(3))
	  vec2(1) = vec2(1) - xpos(i+vec(3)-1)
	  vec2(2) = vec2(2) - ypos(i+vec(3)-1)
	  vec2(3) = vec2(3) - zpos(i+vec(3)-1)
	  mag = sqrt(vec2(1)*vec2(1) + vec2(2)*vec2(2) + vec2(3)*vec2(3))
	  vec2(1) = vec2(1) / mag
	  vec2(2) = vec2(2) / mag
	  vec2(3) = vec2(3) / mag

	  ! Average of vectors vec1 and vec2 give the point which our pointing vector runs through (from O)
	  ! So, pointing vector is |t+p| which is the (H1-H2)->O vector
	  pvec(1) = vec1(1) + vec2(1)
	  pvec(2) = vec1(2) + vec2(2)
	  pvec(3) = vec1(3) + vec2(3)

	  ! Normalise pointing vector
	  mag = sqrt(pvec(1)*pvec(1) + pvec(2)*pvec(2) + pvec(3)*pvec(3))
	  pvec(1) = pvecmult * pvec(1) / mag
	  pvec(2) = pvecmult * pvec(2) / mag
	  pvec(3) = pvecmult * pvec(3) / mag

	  ! Calculate centroid of vector to determine position on/above surface
	  tx = xpos(i+vec(3)-1) - pvec(1) * 0.5
	  ty = ypos(i+vec(3)-1) - pvec(2) * 0.5
	  tz = zpos(i+vec(3)-1) - pvec(3) * 0.5

	  ! Use COM z instead of vector-z if necessary
	  if (usecom) tz = comz(targetsp,n)

	  ! Get x/y position on surface
	  ! Fold into -0.5->0.5
	  if (tx.lt.(-0.5*cell(1))) tx = tx + cell(1)
	  if (tx.gt.(0.5*cell(1))) tx = tx - cell(1)
	  if (ty.lt.(-0.5*cell(5))) ty = ty + cell(5)
	  if (ty.gt.(0.5*cell(5))) ty = ty - cell(5)
	  if (tz.lt.(-0.5*cell(9))) tz = tz + cell(9)
	  if (tz.gt.(0.5*cell(9))) tz = tz - cell(9)
	!write(0,*) patchz, tz
	  ! Check z-distance
	  !write(0,*) dabs(tz-patchz), patchz, minz,maxz
	  if ((dabs(tz - patchz).gt.maxz).or.(dabs(tz-patchz).lt.minz)) then
	    i = i + s_natoms(targetsp)
	    cycle
	  end if
	  !write(34,*) patchz, tz, tz - patchz
	  ! Get bin position on surface and pointing vector bins
	  xbin = tx / xdelta
	  ybin = ty / ydelta
	  zbin = int(pvec(3) / 0.01)
	  ! Get angle between surf-O-H1 and surf-O-H2. Already have the O-H unit vectors so calculate angle made with {0,0,-1}
	  anglebin1 = int(acos(vec1(3)*(-1.0))*radcon)
	  anglebin2 = int(acos(vec2(3)*(-1.0))*radcon)
	  ! Calculate velocity in XY plane (OXYGEN ATOM ONLY)
	  vel = sqrt(xvel(i+vec(3)-1) * xvel(i+vec(3)-1) + yvel(i+vec(3)-1) * yvel(i+vec(3)-1))

	  ! Store data
	  map(xbin,ybin,1) = map(xbin,ybin,1) + pvec(1)
	  map(xbin,ybin,2) = map(xbin,ybin,2) + pvec(2)
	  map(xbin,ybin,3) = map(xbin,ybin,3) + pvec(3)
	  map(xbin,ybin,4) = map(xbin,ybin,4) + (tz - patchz)
	  map(xbin,ybin,5) = map(xbin,ybin,5) + vel
	  histo(zbin) = histo(zbin) + 1.0
	  if (doangles) anglehist(xbin,ybin,1,anglebin1) = anglehist(xbin,ybin,1,anglebin1) + 1
	  if (doangles) anglehist(xbin,ybin,2,anglebin2) = anglehist(xbin,ybin,2,anglebin2) + 1
	  !write(66,*) vec1(3),patchz,dabs(vec1(3) - patchz)
	  acc(xbin,ybin) = acc(xbin,ybin) + 1.0
	  ! Accumulate global average
	  globalsum = globalsum + (vec1(3) - patchz)
	  globalsumpos = globalsumpos + vec1(3)
	  globalacc = globalacc + 1.0

	  i = i + s_natoms(targetsp)
	end do
	
	! Next frame (or finish)
116	if (nframes.EQ.endf) goto 801 
	goto 102

119	write(0,*) "HISTORY file ended prematurely!"
	write(0,"(A,I5,A)") "Managed ",nframesused," frames before error."

	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "Problem with HISTORY file - failed to read header."
	goto 999
800	write(0,*) "End of unformatted HISTORY file found."
801	write(0,*) "Finished."

	! Write out data
	open(unit=11,form="formatted",file=basename(1:baselen)//"map"//CHAR(48+targetsp)//".vf",status="replace")
	open(unit=12,form="formatted",file=basename(1:baselen)//"zmap"//CHAR(48+targetsp)//".surf",status="replace")
	open(unit=13,form="formatted",file=basename(1:baselen)//"zmap"//CHAR(48+targetsp)//".vf",status="replace")
	open(unit=14,form="formatted",file=basename(1:baselen)//"map"//CHAR(48+targetsp)//".surf",status="replace")
	open(unit=15,form="formatted",file=basename(1:baselen)//"map"//CHAR(48+targetsp)//".ang1",status="replace")
	open(unit=16,form="formatted",file=basename(1:baselen)//"map"//CHAR(48+targetsp)//".ang2",status="replace")
	open(unit=17,form="formatted",file=basename(1:baselen)//"vel"//CHAR(48+targetsp)//".surf",status="replace")
	write(12,*) 2*ngrid+1, 2*ngrid+1
	write(12,"(6f8.4)") xdelta,0.0,0.0,0.0,ydelta,0.0
	write(12,"(3f10.4)") -0.5*cell(1),-0.5*cell(5),0.0
	write(12,*) "xyz"
	write(14,*) 2*ngrid+1, 2*ngrid+1
	write(14,"(6f8.4)") xdelta,0.0,0.0,0.0,ydelta,0.0
	write(14,"(3f10.4)") -0.5*cell(1),-0.5*cell(5),0.0
	write(14,*) "xyz"
	write(17,*) 2*ngrid+1, 2*ngrid+1
	write(17,"(6f8.4)") xdelta,0.0,0.0,0.0,ydelta,0.0
	write(17,"(3f10.4)") -0.5*cell(1),-0.5*cell(5),0.0
	write(17,*) "xyz"
	do ybin=-ngrid,ngrid
	  do xbin=-ngrid,ngrid
	    if (acc(xbin,ybin).gt.0.5) then
	      map(xbin,ybin,1) = map(xbin,ybin,1) / acc(xbin,ybin)
	      map(xbin,ybin,2) = map(xbin,ybin,2) / acc(xbin,ybin)
	      map(xbin,ybin,3) = map(xbin,ybin,3) / acc(xbin,ybin)
	      map(xbin,ybin,4) = map(xbin,ybin,4) / acc(xbin,ybin)
	      map(xbin,ybin,5) = map(xbin,ybin,5) / acc(xbin,ybin)
	    else
	      ! Set 'dead' points to the average surface value
	!write(0,*) globalsum,globalacc
              map(xbin,ybin,3) = -2.0
              map(xbin,ybin,4) = globalsum / globalacc
              map(xbin,ybin,5) = -1.0
	    end if
	    ! Average pointing vectors as function over surface
	    write(11,"(6f12.5)") xbin*xdelta, ybin*ydelta, 0.0, (map(xbin,ybin,n),n=1,3)
	    ! Surface zmap relative to calculation centre (patchz)
	    write(12,"(2f12.5)") map(xbin,ybin,4),acc(xbin,ybin)
	    ! Average pointing vectors as function over surface, including height component relative to patch
	    write(13,"(6f12.5)") xbin*xdelta, ybin*ydelta, map(xbin,ybin,4), (map(xbin,ybin,n),n=1,3)
	    ! Map of averages of z-components of vectors over surface
	    write(14,"(f10.5)") map(xbin,ybin,3)
	    ! Histogram of surf-O-H1 angles
	    if (doangles) write(15,"(f10.5)") anglehist(xbin,ybin,1,0:180)
	    ! Histogram of surf-O-H2 angles
	    if (doangles) write(16,"(f10.5)") anglehist(xbin,ybin,2,0:180)
	    ! Velocity average over surface
	    write(17,"(f10.5)") map(xbin,ybin,5)
	  end do
	end do
	close(11)
	close(12)
	close(13)
	close(14)
	close(15)
	close(16)
	close(17)
	! Write histogram data
	open(unit=15,form="formatted",file=basename(1:baselen)//"zhist"//CHAR(48+targetsp)//".dat",status="replace")
	do zbin=-100,100
	  write(15,"(2f12.5)") zbin*0.01, histo(zbin) / globalacc
	end do
	close(15)

	write(0,*) "Finished!"
999	close(10)
	close(20)

	end program point
