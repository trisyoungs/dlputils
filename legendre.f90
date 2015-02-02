!	** legendre.f90 **
!	Calculates the first and second-order Legendre polynomial quantities over x, y, and z-bins

	program legendre
	use dlprw; use utility
	implicit none
	integer :: maxframes,vec(3),baselen,nargs,success,n,i,j,k,targetsp,l
	integer :: nframes, nframesused, bin(3), zbin, nzbins, ngridbins(3)
	integer :: framestodiscard=-1, nskipped
	character*80 :: hisfile,outfile,basename,temp,altheaderfile
	character*8 :: discard
	logical :: altheader = .FALSE., usecom = .FALSE.
	real*8,allocatable :: accxyz(:,:,:,:), legxyz(:,:,:,:), legz(:,:), accz(:,:)
	real*8 :: vec1(3), vec2(3), pvec(3), pos(3), mag, tx, ty, tz, normal(3), binwidth, origin(3), gridspacing(3)
	real*8, parameter :: radcon = 57.29577951d0
	integer, parameter :: nlegs = 2
	real*8 :: leg(nlegs), total
	integer :: iargc

	write(0,*) "*** legendre"

	nargs = iargc()
	if (nargs.LT.2) then
	  write(*,"(a)") "Usage: legendre <HISTORYfile> <OUTPUTfile> [-options]"
	  write(*,"(a)") "        [-target sp]            Set species sp to be the target"
	  write(*,"(a)") "        [-vector i1 i2 j]       Atoms to use for vector calculation in target species: v = | (i1/i2)->j | = +1"
	  write(*,"(a)") "        [-normal x y z]         Define surface normal (default = {0,0,1})"
	  write(*,"(a)") "        [-origin x y z]         Use specified coordinate as origin/MIM point (default = {0.0,0.0,0.0})"
	  write(*,"(a)") "        [-spacing dx dy dx]     Use spacings between gridpoints for grid map (default = {0.5,0.5,0.5})"
	  write(*,"(a)") "        [-binwidth dz]          Use binwidth for z-dependent average (default = 0.1)"
	  write(*,"(a)") "        [-header file]          Use specified file to get header"
	  write(*,"(a)") "        [-discard n]               Skip n frames at start of trajectory"
	  write(*,"(a)") "        [-com]                  Use centre of mass for z-coordinate check rather than vector centre"
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
	vec(1) = 1
	vec(2) = 2
	vec(3) = 3
	binwidth = 0.1
	normal = (/ 0.0, 0.0, 1.0 /)
	gridspacing = 0.5
	origin = 0.0
	targetsp = 1

	n = 2
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-vector") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") vec(1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") vec(2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") vec(3)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Vector for species ",targetsp," calculated from: (",vec(1),",", &
	        & vec(2),")->",vec(3)
	      write(20,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Vector for species ",targetsp," calculated from: (",vec(1),",", &
	        & vec(2),")->",vec(3)
	    case ("-target")
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") targetsp
	      write(0,"(A,I4)") "Target species is = ",targetsp
	      write(20,"(A,I4)") "Target species is = ",targetsp
	    case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
              write(20,"(A)") "Alternative header file supplied."
              altheader = .TRUE.
	    case ("-origin") 
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") origin(1)
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") origin(2)
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") origin(3)
	      write(0,"(A,3f12.6)") "Coordinate origin: ", origin(1), origin(2), origin(3)
	      write(20,"(A,3f12.6)") "Coordinate origin: ", origin(1), origin(2), origin(3)
	    case ("-spacing") 
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") gridspacing(1)
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") gridspacing(2)
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") gridspacing(3)
	      write(0,"(A,3f12.6)") "3D map gridspacing: ", gridspacing(1), gridspacing(2), gridspacing(3)
	      write(20,"(A,3f12.6)") "3D map gridspacing: ", gridspacing(1), gridspacing(2), gridspacing(3)
	    case ("-binwidth")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") binwidth
	      write(0,"(A,f8.4)") "Bin width for z-dependent properties:", binwidth
	      write(20,"(A,f8.4)") "Bin width for z-dependent properties:", binwidth
	    case ("-com")
	      write(0,"(A)") "Molecule COM will be used for z-coordinate."
	      write(20,"(A)") "Molecule COM will be used for z-coordinate."
	      usecom = .TRUE.
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

	! Determine array sizes required
	ngridbins(1) = cell(1) / gridspacing(1)
	ngridbins(2) = cell(5) / gridspacing(2)
	ngridbins(3) = cell(9) / gridspacing(3)
	nzbins = cell(9) / binwidth + 1
	! Allocate arrays
	allocate(accxyz(nlegs,0:ngridbins(1),0:ngridbins(2),0:ngridbins(3)))
	allocate(legxyz(nlegs,0:ngridbins(1),0:ngridbins(2),0:ngridbins(3)))
	allocate(accz(nlegs,0:nzbins))
	allocate(legz(nlegs,0:nzbins))
	legxyz = 0.0
	accxyz = 0.0
	legz = 0.0
	accz = 0.0

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes = 0
	nskipped = 0
	nframesused = 0

102	success=readframe()
	IF (success.eq.1) goto 801  ! End of file encountered....
	IF (success.eq.-1) goto 119  ! File error....
	nskipped = nskipped + 1
	if (nskipped.le.framestodiscard) goto 102
	nframes=nframes+1
	if (MOD(nframes,100).eq.0) write(0,*) nframes

	nframesused = nframesused + 1

	! Calculate COM if needed
	if (usecom) call calc_com()

	i = s_start(targetsp)
	do n=1,s_nmols(targetsp)

	  ! Get mim unit vector position1 from first atom with third atom (in vec1)
	  call pbc(xpos(i+vec(1)-1),ypos(i+vec(1)-1),zpos(i+vec(1)-1), &
		& xpos(i+vec(3)-1),ypos(i+vec(3)-1),zpos(i+vec(3)-1), tx, ty, tz)
	  pos(1) = tx
	  pos(2) = ty
	  pos(3) = tz
	  vec1(1) = xpos(i+vec(3)-1) - tx
	  vec1(2) = ypos(i+vec(3)-1) - ty
	  vec1(3) = zpos(i+vec(3)-1) - tz
	  mag = sqrt(vec1(1)*vec1(1) + vec1(2)*vec1(2) + vec1(3)*vec1(3))
	  vec1(1) = vec1(1) / mag
	  vec1(2) = vec1(2) / mag
	  vec1(3) = vec1(3) / mag
	  
	  ! Get MIM unit vector position from second atom with third atom (in vec2)
	  call pbc(xpos(i+vec(2)-1),ypos(i+vec(2)-1),zpos(i+vec(2)-1), &
		& xpos(i+vec(3)-1),ypos(i+vec(3)-1),zpos(i+vec(3)-1), tx, ty, tz)
	  pos(1) = pos(1) + tx
	  pos(2) = pos(2) + ty
	  pos(3) = pos(3) + tz
	  vec2(1) = xpos(i+vec(3)-1) - tx
	  vec2(2) = ypos(i+vec(3)-1) - ty
	  vec2(3) = zpos(i+vec(3)-1) - tz
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
	  pvec(1) = pvec(1) / mag
	  pvec(2) = pvec(2) / mag
	  pvec(3) = pvec(3) / mag

	  ! Calculate centroid of vector to determine position on/above surface
	  ! Use COM instead of vector-z if necessary
	  if (usecom) then
	    pos(1) = comx(targetsp,n)
	    pos(2) = comy(targetsp,n)
	    pos(3) = comz(targetsp,n)
	  else
	    pos(1) = (xpos(i+vec(3)-1) + pos(1)) / 3.0
	    pos(2) = (ypos(i+vec(3)-1) + pos(2)) / 3.0
	    pos(3) = (zpos(i+vec(3)-1) + pos(3)) / 3.0
	  end if

	  ! Get position relative to specified origin
	  ! Acquire coordinate relative to origin and with all coords in (fractional) range 0-1
	  pos = pos - origin
	  ! Fold into 0->1
	  call unitfold(pos(1),pos(2),pos(3))
	  !if (pos(1).lt.0.0) pos(1) = pos(1) + int(abs(pos(1)/cell(1))+1)*cell(1)
	  !if (pos(1).gt.cell(1)) pos(1) = pos(1) - cell(1)
	  !if (pos(2).lt.0.0) pos(2) = pos(2) + int(abs(pos(2)/cell(5))+1)*cell(1)
	  !if (pos(2).gt.cell(5)) pos(2) = pos(2) - cell(5)
	  !if (pos(3).lt.0.0) pos(3) = pos(3) + cell(9)
	  !if (pos(3).gt.cell(9)) pos(3) = pos(3) - cell(9)

	  ! Get gridbins and zbin
	  bin = pos / gridspacing
	  zbin = pos(3) / binwidth

	  ! Calculate Legendre polynomials
	  ! 1) First order - just project pointing vector onto surface normal ( = cos(theta))
	  leg(1) = pvec(1)*normal(1) + pvec(2)*normal(2) + pvec(3)*normal(3)
	  ! 2) Second order
	  leg(2) = 0.5 * (3.0 * leg(1)*leg(1) - 1.0)

	  ! Store data
	  do l=1,nlegs
	    legxyz(l,bin(1),bin(2),bin(3)) = legxyz(l,bin(1),bin(2),bin(3)) + leg(l)
	    accxyz(l,bin(1),bin(2),bin(3)) = accxyz(l,bin(1),bin(2),bin(3)) + 1.0
	    legz(l,zbin) = legz(l,zbin) + leg(l)
	    accz(l,zbin) = accz(l,zbin) + 1.0
	  end do

	  ! Next molecule
	  i = i + s_natoms(targetsp)
	end do
	
	! Next frame (or finish)
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
	do l=1,nlegs
	  open(unit=11,form="formatted",file=basename(1:baselen)//"sp"//CHAR(48+targetsp)//".leg"//CHAR(48+l),status="replace")
	  open(unit=12,form="formatted",file=basename(1:baselen)//"sp"//CHAR(48+targetsp)//".leg"//CHAR(48+l)//".pdens",status="replace")
	  if (l.eq.1) open(unit=13,form="formatted",file=basename(1:baselen)//"sp"//CHAR(48+targetsp)//".acc.pdens",status="replace")
	  if (l.eq.1) open(unit=15,form="formatted",file=basename(1:baselen)//"sp"//CHAR(48+targetsp)//".normacc.pdens",status="replace")
	  open(unit=14,form="formatted",file=basename(1:baselen)//"sp"//CHAR(48+targetsp)//".normleg"//CHAR(48+l)//".pdens",status="replace")
	  write(12,*) ngridbins(1)+1, ngridbins(2)+1, ngridbins(3)+1
	  write(12,"(9f6.2)") gridspacing(1),0.0,0.0,0.0,gridspacing(2),0.0,0.0,0.0,gridspacing(3)
	  write(12,"(3f10.4)") 0.0,0.0,0.0
	  write(12,*) "zyx"
	  if (l.eq.1) then
	    write(13,*) ngridbins(1)+1, ngridbins(2)+1, ngridbins(3)+1
	    write(13,"(9f6.2)") gridspacing(1),0.0,0.0,0.0,gridspacing(2),0.0,0.0,0.0,gridspacing(3)
	    write(13,"(3f10.4)") 0.0,0.0,0.0
	    write(13,*) "zyx"
	  end if
	  write(14,*) ngridbins(1)+1, ngridbins(2)+1, ngridbins(3)+1
	  write(14,"(9f6.2)") gridspacing(1),0.0,0.0,0.0,gridspacing(2),0.0,0.0,0.0,gridspacing(3)
	  write(14,"(3f10.4)") 0.0,0.0,0.0
	  write(14,*) "zyx"

	  ! Write density map
	  do i=0,ngridbins(1)
	    do j=0,ngridbins(2)
	      do k=0,ngridbins(3)
	        if (accxyz(l,i,j,k).gt.0.5) then
	          legxyz(l,i,j,k) = legxyz(l,i,j,k) / accxyz(l,i,j,k)
	        end if
	        write(12,"(f12.5)") legxyz(l,i,j,k)
	        if (l.eq.1) write(13,"(f12.5)") accxyz(l,i,j,k) / nframes
	        write(14,"(f12.5)") legxyz(l,i,j,k) * (accxyz(l,i,j,k) / nframes)
	      end do
	    end do
	  end do
	      
	  ! Z-profiles
	  write(11,"(a)") "# Distance    P        Nperframe   SUM(N)"
	  total = 0.0
	  do i=0,nzbins
	    if (accz(l,i).gt.0.5) legz(l,i) = legz(l,i) / accz(l,i)
	    total = total + accz(l,i) / nframesused
	    write(11,"(4f12.5)") i*binwidth, legz(l,i), accz(l,i) / nframesused, total
	  end do

	  ! Close files ready for next Leg
	  close(11)
	  close(12)
	  if (l.eq.1) close(13)
	  close(14)
	end do

	write(0,*) "Finished!"
999	close(10)
	close(20)

	end program legendre
