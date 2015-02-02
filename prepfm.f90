!	** prepfm **
!	Prepare input for Siesta calculations to get forced for force-matching
!	For individual ions, outputs the MIM'd coordinates centred at 0,0,0

	program prepfm
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,outfile
	character*20 :: temp, sysname
	character*6 :: typenames(20)
	character*4 :: frameid
	character*2 :: ionid
	integer :: n,m,baselen,nframes,success,nargs,calctype,offset,frameskip
	integer :: framestodo = -1, dot, space, ntypes, typeelements(20), sysatoms
	integer, allocatable :: s_elements(:)
	real*8 :: mixweight, cog(3), tpos(3)
	integer :: meshcut, eshift
	integer :: iargc
	logical :: found

	nargs = iargc()
	if (nargs.ne.8) then
	  write(0,"(A)") "Usage : prepfm <HISTORYfile> <OUTfile> <out> <meshcut> <eshift> <mixweight> <sysname> <frameskip>"
	  write(0,"(A)") "	 <out> is 'cell' for original config, 'cations' or 'anions'"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp)
	calctype = -1
	select case (temp)
	  case ("cell"); calctype = 0
	  case ("cations"); calctype = 1
	  case ("anions"); calctype = 2
	end select
	if (calctype.eq.-1) stop "Unrecognised calculation type."
	call getarg(4,temp); read(temp,"(i6)") meshcut
	call getarg(5,temp); read(temp,"(i6)") eshift
	call getarg(6,temp); read(temp,"(f10.5)") mixweight
	call getarg(7,sysname)
	call getarg(8,temp); read(temp,"(i6)") frameskip

	! Make up output basename
	dot = 0; space = 0
	do n=1,80
	  if (hisfile(n:n).eq.".") dot = n       ! Pos of last dot
	  if ((hisfile(n:n).eq." ").and.(space.eq.0)) space = n  ! Pos of first space
	 end do
	if (dot.EQ.0) dot=space

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
	! Now, read in the history header so that we have cell()
	if (readheader().EQ.-1) goto 799

	! Convert atomnames to elements (super-crude way...) and Siesta-ise the atomnames
	allocate(s_elements(natms))
	do n=1,natms
	  select case (atmname(n)(1:1))
	    case("H"); s_elements(n) = 1; atmname(n) = "H-PBE"
	    case("B"); s_elements(n) = 5; atmname(n) = "B-PBE"
	    case("C")
	      if (atmname(n)(2:2).eq."l") then
		s_elements(n) = 17; atmname(n) = "Cl-PBE"
	      else
		s_elements(n) = 6; atmname(n) = "C-PBE"
	      end if
	    case("N"); s_elements(n) = 7; atmname(n) = "N-PBE"
	    case("O"); s_elements(n) = 8; atmname(n) = "O-PBE"
	    case("F"); s_elements(n) = 9; atmname(n) = "F-PBE"
	    case("P"); s_elements(n) = 15; atmname(n) = "P-PBE"
	    case("S"); s_elements(n) = 16; atmname(n) = "S-PBE"
	  end select
	end do

	! Construct unique type list
	ntypes = 0
	! Cations... (assumed to be species 1)
	if (calctype.lt.2) then
	  do n=s_start(1),s_natoms(1)
	    ! Check for this element in the list
	    found = .FALSE.
	    do m=1,ntypes
	      if (typeelements(m).eq.s_elements(n)) found = .TRUE.
	    end do
	    if (.not.found) then
	      ntypes = ntypes + 1
	      typeelements(ntypes) = s_elements(n)
	      typenames(ntypes) = atmname(n)
	    end if
	  end do
	end if
	! Anions... (assumed to be species 2)
	if (calctype.ne.1) then
	  do n=s_start(2),s_start(2)+s_natoms(2)
	    ! Check for this element in the list
	    found = .FALSE.
	    do m=1,ntypes
	      if (typeelements(m).eq.s_elements(n)) found = .TRUE.
	    end do
	    if (.not.found) then
	      ntypes = ntypes + 1
	      typeelements(ntypes) = s_elements(n)
	      typenames(ntypes) = atmname(n)
	    end if
	  end do
	end if
	write(0,"(a,i1,a)") "There are ",ntypes," unique types in the system."

	! Convert atom elements held in s_elements into type identifiers
	! Do search by typenames...
	do n=1,natms
	  do m=1,ntypes
	    if (atmname(n).eq.typenames(m)) s_elements(n) = m
	  end do
	end do
	   

100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	if (mod(nframes-1,frameskip).ne.0) goto 101

	write(frameid,"(i4)") nframes
	do n=1,3
	  if (frameid(n:n).eq." ") frameid(n:n) = "0"
	end do

	! Begin output of commands to file
	if (calctype.eq.0) then
	  ! Output entire cells at each step
	  outfile = hisfile(1:dot-1)//"-frame"//frameid//"-cell.fdf"
	  open(unit=15,file=outfile,form='formatted',status='replace')
	  write(15,"(a,a30)") "SystemName          ",dlpheader
	  write(15,"(a,a30)") "SystemLabel         ",sysname
	  write(15,"(a,i4)")  "NumberOfAtoms       ",natms
	  write(15,"(a,i4)")  "NumberOfSpecies     ",ntypes
	  call writetypes(ntypes,typenames,typeelements)
	  call writebulk(calctype,meshcut,eshift,mixweight,cell(1))
	  ! Write atomic coordinates
	  ! s_elements now holds the atom types
	  write(15,"(a)") "AtomicCoordinatesFormat  Ang"
	  write(15,"(a)") "%block AtomicCoordinatesAndAtomicSpecies"
	  do n=1,natms
	    write(15,"(3(f12.8),2x,i2)") xpos(n),ypos(n),zpos(n),s_elements(n)
	  end do
	  write(15,"(a)") "%endblock AtomicCoordinatesAndAtomicSpecies"
	  close(15)
	else if (calctype.eq.1) then
	  ! Output cations at each step
	  offset = s_start(1)
	  do n=1,s_nmols(1)
	    write(ionid,"(i2)") n
	    if (ionid(1:1).eq." ") ionid(1:1) = "0"
	    outfile = hisfile(1:dot-1)//"-frame"//frameid//"-cation"//ionid//".fdf"
	    open(unit=15,file=outfile,form='formatted',status='replace')
	    write(15,"(a,a30)") "SystemName          ",dlpheader
	    write(15,"(a,a30)") "SystemLabel         ",sysname
	    write(15,"(a,i4)")  "NumberOfAtoms       ",s_natoms(1)
	    write(15,"(a,i4)")  "NumberOfSpecies     ",ntypes
	    call writetypes(ntypes,typenames,typeelements)
	    call writebulk(calctype,meshcut,eshift,mixweight,0.0)
	    ! Write atomic coordinates
	    ! s_elements now holds the atom types
	    write(15,"(a)") "AtomicCoordinatesFormat  Ang"
	    write(15,"(a)") "%block AtomicCoordinatesAndAtomicSpecies"
	    ! Must get minimim image coordinates, centred at zero.
	    ! So, for this cation calculate mim coords w.r.t. first atom, calculating COG at same time
	    cog(1) = xpos(offset)
	    cog(2) = ypos(offset)
	    cog(3) = zpos(offset)
	    do m=offset+1,offset+s_natoms(1)-1
	      call pbc(xpos(m),ypos(m),zpos(m),xpos(offset),ypos(offset),zpos(offset),tpos(1),tpos(2),tpos(3))
	      cog = cog + tpos
	      xpos(m) = tpos(1)
	      ypos(m) = tpos(2)
	      zpos(m) = tpos(3)
	    end do
	    cog = cog / s_natoms(1)
	    do m=offset,offset+s_natoms(1)-1
	      write(15,"(3(f12.8),2x,i2)") xpos(m)-cog(1),ypos(m)-cog(2),zpos(m)-cog(3),s_elements(m)
	    end do
	    write(15,"(a)") "%endblock AtomicCoordinatesAndAtomicSpecies"
	    close(15)
	    offset = offset + s_natoms(1)
	  end do
	else if (calctype.eq.2) then
	  ! Output anions at each step
	  offset = s_start(2)
	  do n=1,s_nmols(2)
	    write(ionid,"(i2)") n
	    if (ionid(1:1).eq." ") ionid(1:1) = "0"
	    outfile = hisfile(1:dot-1)//"-frame"//frameid//"-anion"//ionid//".fdf"
	    open(unit=15,file=outfile,form='formatted',status='replace')
	    write(15,"(a,a30)") "SystemName          ",dlpheader
	    write(15,"(a,a30)") "SystemLabel         ",sysname
	    write(15,"(a,i4)")  "NumberOfAtoms       ",s_natoms(2)
	    write(15,"(a,i4)")  "NumberOfSpecies     ",ntypes
	    call writetypes(ntypes,typenames,typeelements)
	    call writebulk(calctype,meshcut,eshift,mixweight,0.0)
	    ! Write atomic coordinates
	    ! s_elements now holds the atom types
	    write(15,"(a)") "AtomicCoordinatesFormat  Ang"
	    write(15,"(a)") "%block AtomicCoordinatesAndAtomicSpecies"
	    ! Must get minimim image coordinates, centred at zero.
	    ! So, for this cation calculate mim coords w.r.t. first atom, calculating COG at same time
	    cog(1) = xpos(offset)
	    cog(2) = ypos(offset)
	    cog(3) = zpos(offset)
	    do m=offset+1,offset+s_natoms(2)-1
	      call pbc(xpos(m),ypos(m),zpos(m),xpos(offset),ypos(offset),zpos(offset),tpos(1),tpos(2),tpos(3))
	      cog = cog + tpos
	      xpos(m) = tpos(1)
	      ypos(m) = tpos(2)
	      zpos(m) = tpos(3)
	    end do
	    cog = cog / s_natoms(1)
	    do m=offset,offset+s_natoms(2)-1
	      write(15,"(3(f12.8),2x,i2)") xpos(m)-cog(1),ypos(m)-cog(2),zpos(m)-cog(3),s_elements(m)
	    end do
	    write(15,"(a)") "%endblock AtomicCoordinatesAndAtomicSpecies"
	    close(15)
	    offset = offset + s_natoms(2)
	  end do
	end if

	if (nframes.EQ.framestodo) goto 800
	! Next frame
	goto 101

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended before framestodo was fulfilled..."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",framestodo,")"
	goto 801
800	write(0,*) "Framestodo was fulfilled."
801	write(0,*) ""

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program prepfm

	subroutine writetypes(ntypes,names,elements)
	implicit none
	character*6 :: names(20)
	integer :: elements(20),ntypes,n
	write(15,*) ""
	write(15,"(a)")  "%block ChemicalSpeciesLabel"
	do n=1,ntypes
	  write(15,"(i2,i3,2x,a8)") n,elements(n),names(n)
	end do
	write(15,"(a)")  "%endblock ChemicalSpeciesLabel"
	write(15,*) ""
	end subroutine writetypes

	subroutine writebulk(calctype,meshcut,eshift,mixweight,lattice)
	implicit none
	integer :: meshcut, eshift,calctype
	real*8 :: mixweight, lattice
	write(15,"(a,i3,a)") "PAO.EnergyShift     ",eshift," meV"
	write(15,"(a,i3,a)") "MeshCutoff          ",meshcut," Ry"
	write(15,"(a,f4.1,a)") "DM.MixingWeight     ",mixweight
	if (calctype.eq.0) write(15,"(a,f14.10,a)") "LatticeConstant     ",lattice," Ang"
	if (calctype.eq.1) write(15,"(a)") "NetCharge           1"
	if (calctype.eq.2) write(15,"(a)") "NetCharge           -1"
	write(15,*) ""
	write(15,"(a)") "PAO.BasisSize       DZP"
	write(15,"(a)") "XC.functional       GGA"
	write(15,"(a)") "XC.authors          PBE"
	write(15,"(a)") "DM.NumberPulay      4"
	write(15,"(a)") "DM.Tolerance        1.d-6"
	write(15,"(a)") "SolutionMethod      diagon"
	write(15,"(a)") "MaxSCFIterations    100"
	write(15,"(a)") "ParallelOverK       true"
	write(15,"(a)") "MD.TypeOfRun        Verlet"
	write(15,"(a)") "MD.InitialTemperature  0.0  K"
	write(15,"(a)") "MD.InitialTimeStep    1"
	write(15,"(a)") "MD.FinalTimeStep    1"
	write(15,"(a)") "WriteCoorXmol       true"
	write(15,"(a)") "WriteForces         true"
	write(15,*) ""
	end subroutine writebulk
