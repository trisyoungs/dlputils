	! Routines to read/write DL_POLY formatted and unformatted history files.
	! Also, routine to probe a DL_POLY OUTPUT file for molecule/species information.

	module dlprw
	  ! Limits
	  integer, parameter :: MAXSP = 10
	  ! History file
	  character*80 :: dlpfn_out, dlpfn_his
	  integer :: his_in_format = -1, dlpun_out, dlpun_his, dlpun_cfg, dlpun_newhis, his_out_format = -1
	  ! History file data
	  character*80 :: dlpheader
	  real*8,allocatable,dimension(:) :: xpos,ypos,zpos,charge,mass,xvel,yvel,zvel, &
		& xfor,yfor,zfor
	  real*8 :: cell(9), tstep
	  logical :: has_icell = .false.
	  character*8, allocatable :: atmname(:)
	  integer :: keytrj, imcon, natms, nstep
	  ! Outfile data (species composition, count etc)
          character*20, allocatable :: s_name(:)
          integer :: nspecies,maxmols,maxatoms
	  integer, allocatable :: s_natoms(:),s_nmols(:),s_start(:),s_molstart(:)
	  integer :: s_totalMols
	  ! Config file data
	  real*8, allocatable,dimension(:,:) :: cfgxyz,cfgvel,cfgfrc
	  real*8 :: cfgcell(9)
	  character*8, allocatable :: cfgname(:)
	  character*80 :: cfgheader
	  ! Flag
	  logical :: xyz_alloc = .false.
	contains

!	\/\/\/\/\/\/\/
!	** openhis **
!	/\/\/\/\/\/\/\
!	Determine the format of the history file (if possible, otherwise ask).

	subroutine openhis(fname,uno)
	implicit none
	character*80 :: fname,tempc
	integer :: uno,tempi,baselen
	dlpfn_his = fname; dlpun_his = uno
        ! Determine the FORM of the HISTORY file, in an horrifically crude way because
	! the INQUIRE statement is a bit 'chocolate fireguard'.
	! Check last four characters of filename....
	! Ascertain length of basename....
	do baselen=80,1,-1
	  if (fname(baselen:baselen).ne." ") exit
	end do
	if (baselen.gt.4) then
	  if (fname(baselen-3:baselen).eq."HISf") goto 5
	  if (fname(baselen-3:baselen).eq."HISu") goto 10
	endif
	! First, open the file as FORMATTED and attempt to read the first integers...
	open(file=dlpfn_his,unit=40,form='formatted',err=10,status='old')
	read(40,"(A80)",err=10) tempc
	read(40,"(2I10)",err=10) tempi,tempi
	! If we get to here then it *is* a formatted file.
5	his_in_format = 1; close(40)
	open(unit=dlpun_his,file=dlpfn_his,form='formatted',status='old')
	write(0,*) "History file is FORMATTED."
	return
	! Howeverm if we arrive here then the read failed - thus its UNformatted (hopefully)
10	his_in_format = 0; close(40)
	open(unit=dlpun_his,file=dlpfn_his,form='unformatted',status='old')
	write(0,*) "History file is UNFORMATTED."
	return
	write(0,*) "The specified history file does not exist!!"
	write(0,"(A,A)") "--> ",dlpfn_his
	stop
	end subroutine openhis


!	\/\/\/\/\/\/\/\/
!	** readframe **
!	/\/\/\/\/\/\/\/\
!	Subroutines to read in a DL_POLY history frame.
!	 --> readframe() and readheader() routines.

	integer function readframe()
	use parse
	implicit none
!	dlpun_his = Fortran unit number to read from.
!	his_in_format = Form of the history file: 0 = Unformatted, 1 = Formatted.
	character*8 :: discard
	integer :: n, discardn
	real*8 :: x1,x2,x3,x4
!12	FORMAT (3g12.4)
!12	FORMAT (g16.4,4x,g16.4,4x,g16.4)    ! For DL_POLY3
14	FORMAT (A8,I10,2F12.6)
	if (his_in_format.eq.-1) stop "HISfile format not yet determined / file not yet opened!"	
	if (his_in_format.eq.0) then   ! Read unformatted history file.....
	  read(dlpun_his,end=91,ERR=92) x1,x2,x3,x4,tstep
	  nstep=IDNINT(x1)
	  keytrj=IDNINT(x3)
	  imcon=IDNINT(x4)
	  if (imcon.GT.0) read(dlpun_his,ERR=93,end=93) cell
	  read(dlpun_his,end=94,ERR=94) (xpos(n),n=1,natms)
	  read(dlpun_his,end=95,ERR=95) (ypos(n),n=1,natms)
	  read(dlpun_his,end=96,ERR=96) (zpos(n),n=1,natms)
	  if (keytrj.GT.0) then
	    read(dlpun_his,end=97,ERR=97) (xvel(n),n=1,natms)
	    read(dlpun_his,end=98,ERR=98) (yvel(n),n=1,natms)
	    read(dlpun_his,end=99,ERR=99) (zvel(n),n=1,natms)
	  end if
	  if (keytrj.GT.1) then
	    read(dlpun_his,end=100,ERR=100) (xfor(n),n=1,natms)
	    read(dlpun_his,end=101,ERR=101) (yfor(n),n=1,natms)
	    read(dlpun_his,end=102,ERR=102) (zfor(n),n=1,natms)
	  end if
	else     ! Read formatted history file....
          !read(dlpun_his,"(8X,4I10,F12.6)",end=91,ERR=103) nstep, natms, keytrj, imcon, tstep
	  if (.not.readline(dlpun_his)) goto 105
	  nstep = argi(2)
	  natms = argi(3)
	  keytrj = argi(4)
	  imcon = argi(5)
	  !write(0,*) nstep, natms, keytrj, imcon, tstep
          if (imcon.GT.0) then     ! Read in the cell coords if necessary
	    !write(0,*) imcon
            !read(dlpun_his,12,ERR=104,end=104) cell
	    !write(0,"(3f12.4)") cell
	    !read(dlpun_his,"(3x,3f14.4)",ERR=104,end=104) cell	! For Agilio's HISTORY files
	    if (.not.readline(dlpun_his)) goto 104
	    cell(1) = argr(1); cell(2) = argr(2); cell(3) = argr(3)
	    if (.not.readline(dlpun_his)) goto 104
	    cell(4) = argr(1); cell(5) = argr(2); cell(6) = argr(3)
	    if (.not.readline(dlpun_his)) goto 104
	    cell(7) = argr(1); cell(8) = argr(2); cell(9) = argr(3)
          end if
          ! Atomic positions/velocities/forces....
          do n=1,natms
	    ! Atom name, mass, charge etc
	    if (.not.readline(dlpun_his)) goto 105
	    atmname(n) = arg(1)(1:8); mass(n) = argr(3); charge(n) = argr(4);
	    if (.not.readline(dlpun_his)) goto 106
	    xpos(n) = argr(1); ypos(n) = argr(2); zpos(n) = argr(3)
	    if (keytrj.gt.0) then
	      if (.not.readline(dlpun_his)) goto 107
	      xvel(n) = argr(1); yvel(n) = argr(2); zvel(n) = argr(3)
	    end if
	    if (keytrj.gt.1) then
	      if (.not.readline(dlpun_his)) goto 108
	      xfor(n) = argr(1); yfor(n) = argr(2); zfor(n) = argr(3)
	    end if
            !read(dlpun_his,14,ERR=105,end=105) discard,discardn,mass(n),charge(n)
            !read(dlpun_his,12,ERR=106,end=106) xpos(n),ypos(n),zpos(n)
	    !if (keytrj.GT.0) read(dlpun_his,12,ERR=107,end=107) xvel(n),yvel(n),zvel(n)
	    !if (keytrj.GT.1) read(dlpun_his,12,ERR=108,end=108) xfor(n),yfor(n),zfor(n)
          end do
	end if

	readframe=0  ; return	! Success
	! Other return conditions:
91	readframe=1  ; return	! No error, but end of file encountered.
92	readframe=-1 ; return	! U - frame data
93	readframe=-2 ; return	! U - cell
94	readframe=-3 ; return	! U - xpositions
95	readframe=-4 ; return	! U - ypositions
96	readframe=-5 ; return	! U - zpositions
97	readframe=-6 ; return	! U - xvelocities
98	readframe=-7 ; return	! U - yvelocities
99	readframe=-8 ; return	! U - zvelocities
100	readframe=-9 ; return	! U - xforces
101	readframe=-10; return	! U - yforces
102	readframe=-11; return	! U - zforces
103	readframe=-12; return	! F - frame data
104	readframe=-13; return	! F - cell
105	readframe=-14; return	! F - atomname,mass,charge
106	readframe=-15; return	! F - positions
107	readframe=-16; return	! F - velocities
108	readframe=-17; return	! F - forces
	end function readframe

	integer function readheader()
	use parse
	implicit none
!	his_in_format = Form of the history file: 0 = Unformatted, 1 = Formatted.
	character*8 discard
	integer :: n, discardn
	real*8 :: discardr,x1,x2,x3,x4
10	FORMAT (A80)
!12	FORMAT (3g12.4)
12	FORMAT (g16.4,4x,g16.4,4x,g16.4)    ! For DL_POLY3
14	FORMAT (A8,I10,2F12.6)
	if (his_in_format.EQ.-1) stop "HISfile format not yet determined / file not yet opened!"
	if (his_in_format.EQ.0) then   ! Read unformatted history file.....
	  write(0,*) "Reading header from unformatted history file..."
	  ! File dlpheader
	  read(dlpun_his,ERR=91,end=91) dlpheader
	  ! Number of atoms, followed by individual names, masses and charges
	  read(dlpun_his,ERR=92,end=92) x1
	  natms=IDNINT(x1)
	  call alloc_xyz
	  read(dlpun_his,ERR=93,end=93) (atmname(n),n=1,natms)
	  read(dlpun_his,ERR=94,end=94) (mass(n),n=1,natms)
	  read(dlpun_his,ERR=95,end=95) (charge(n),n=1,natms)
	  ! Need to spin on and read ahead to get cell(), then come back....
	  read(dlpun_his,end=96,ERR=96) x1,x2,x3,x4,tstep
	  nstep=IDNINT(x1)
	  keytrj=IDNINT(x3)
	  imcon=IDNINT(x4)
	  if (imcon.GT.0) read(dlpun_his,ERR=97,end=97) cell
	  ! Now rewind...
	  REWIND(dlpun_his)
	  read(dlpun_his,ERR=98,end=98) dlpheader
	  read(dlpun_his,ERR=99,end=99) x1
	  read(dlpun_his,ERR=100,end=100) (atmname(n),n=1,natms)
	  read(dlpun_his,ERR=101,end=101) (mass(n),n=1,natms)
	  read(dlpun_his,ERR=102,end=102) (charge(n),n=1,natms)
	else     ! Read formatted history file....
	  write(0,*) "Reading header from formatted history file..."
	  ! File header
	  read(dlpun_his,10,ERR=103,end=103) dlpheader
	  if (dlpheader(1:8).eq."timestep") then
	    write(0,*) "No HISf header records found, but 'timestep' located on first line..."
	    rewind(dlpun_his)
	  else
	    read(dlpun_his,"(2I10)",ERR=104,end=104) keytrj,imcon
	  end if
	  ! Cell params
	  !read(dlpun_his,"(8X,4I10,F12.6)",ERR=105,end=105) nstep, natms, keytrj, imcon, tstep
	  if (.not.readline(dlpun_his)) goto 105
	  nstep = argi(2)
	  natms = argi(3)
	  keytrj = argi(4)
	  imcon = argi(5)
	  ! Read in the cell coords if necessary
	  if (imcon.GT.0) then
	    !read(dlpun_his,12,ERR=106,end=106) cell
	    !read(dlpun_his,"(3x,3f14.4)",ERR=106,end=106) cell	! For Agilio's HISTORY files
	    if (.not.readline(dlpun_his)) goto 106
	    cell(1) = argr(1); cell(2) = argr(2); cell(3) = argr(3)
	    if (.not.readline(dlpun_his)) goto 106
	    cell(4) = argr(1); cell(5) = argr(2); cell(6) = argr(3)
	    if (.not.readline(dlpun_his)) goto 106
	    cell(7) = argr(1); cell(8) = argr(2); cell(9) = argr(3)
          end if
          ! Atomic positions/velocities/forces....
	  call alloc_xyz
          do n=1,natms
	    ! Atom name, mass, charge etc
	    if (.not.readline(dlpun_his)) goto 107
	    atmname(n) = arg(1)(1:8); mass(n) = argr(3); charge(n) = argr(4);
	    if (.not.readline(dlpun_his)) goto 108
	    if (keytrj.gt.0) then
	      if (.not.readline(dlpun_his)) goto 109
	    end if
	    if (keytrj.gt.1) then
	      if (.not.readline(dlpun_his)) goto 110
	    end if
	    !read(dlpun_his,14,ERR=107,end=107) atmname(n),discardn,mass(n),charge(n)
	    !read(dlpun_his,12,ERR=108,end=108) discardr,discardr,discardr
	    !if (keytrj.GT.0) read(dlpun_his,12,ERR=109,end=109) discardr,discardr,discardr
	    !if (keytrj.GT.1) read(dlpun_his,12,ERR=110,end=110) discardr,discardr,discardr
	  end do
	  ! Now, put the hisfile back to the start of the frames data....
	  rewind(dlpun_his)
	  read(dlpun_his,10) dlpheader
	  if (dlpheader(1:8).eq."timestep") then
	    rewind(dlpun_his)
	  else
	    read(dlpun_his,"(2I10)",ERR=104,end=104) keytrj,imcon
	  end if
	  !write(0,*) "At end of header read, last discarded line was:"
	  !write(0,*) discard
	end if
	! Return conditions:
	readheader=0  ; return ! Success
91	readheader=-1 ; return  ! U - header string
92	readheader=-2 ; return	! U - natoms
93	readheader=-3 ; return	! U - atmname
94	readheader=-4 ; return	! U - mass
95	readheader=-5 ; return	! U - charge
96	readheader=-6 ; return	! U - frame - tstep, natoms etc.
97	readheader=-7 ; return	! U - frame - cell
98	readheader=-8 ; return	! U - header (2)
99	readheader=-9 ; return	! U - natoms (2)
100	readheader=-10; return	! U - atmname (2)
101	readheader=-11; return	! U - mass (2)
102	readheader=-12; return	! U - charge (2)
103	readheader=-13; return	! F - header
104	readheader=-14; return	! F - keytrj,imcon
105	readheader=-15; return	! F - nstep,natoms,keytrj,imcon,tstep
106	readheader=-16; return	! F - frame - cell
107	readheader=-17; return	! F - frame - atomname,mass,charge
108	readheader=-18; return	! F - frame - positions
109	readheader=-19; return	! F - frame - velocities
110	readheader=-20; return	! F - frame - forces
	return
	end function readheader

!	\/\/\/\/\/\/\/\/
!	** writeframe **
!	/\/\/\/\/\/\/\/\
!	Subroutines to write out a DL_POLY history frame.
!	--> writeframe and writeheader functions.

	integer function writeframe()
	implicit none
!	dlpun_his = Fortran unit number to read from.
!	his_out_format = Format of history file: 0 = Unformatted, 1 = Formatted.
	integer :: n
12	FORMAT (3e20.12)	 ! coordinate, velocity, or force line
13	FORMAT (A8,4I10,2F12.6)	 ! 'timestep' line
14	FORMAT (A8,I10,2F12.6)	 ! 'atom' line
15	FORMAT (3g20.12)	 ! 'cell' lines
	if (his_out_format.EQ.-1) stop "HISfile format not yet determined / file not yet opened!"
	if (his_out_format.EQ.0) then   ! Write unformatted history file.....
	  write(dlpun_newhis,ERR=91) DBLE(nstep),DBLE(NATMS),DBLE(keytrj),DBLE(imcon),tstep
	  if (imcon.GT.0) write(dlpun_newhis,ERR=91) cell
	  ! Atomic positions/velocities/forces....
	  write(dlpun_newhis,ERR=91) (xpos(n),n=1,natms)
	  write(dlpun_newhis,ERR=91) (ypos(n),n=1,natms)
	  write(dlpun_newhis,ERR=91) (zpos(n),n=1,natms)
	  if (keytrj.GT.0) then
	    write(dlpun_newhis,ERR=91) (xvel(n),n=1,natms)
	    write(dlpun_newhis,ERR=91) (yvel(n),n=1,natms)
	    write(dlpun_newhis,ERR=91) (zvel(n),n=1,natms)	
	  end if
	  if (keytrj.GT.1) then
	    write(dlpun_newhis,ERR=91) (xfor(n),n=1,natms)
	    write(dlpun_newhis,ERR=91) (yfor(n),n=1,natms)
	    write(dlpun_newhis,ERR=91) (zfor(n),n=1,natms)	
	  end if
	else     ! Write formatted history file....
	  write(dlpun_newhis,13,ERR=91) "timestep",nstep,natms,keytrj,imcon,tstep
	  if (imcon.GT.0) write(dlpun_newhis,15,ERR=91) cell
	  do n=1,natms
	    write(dlpun_newhis,14,ERR=91) atmname(n),n,mass(n),charge(n)
	    write(dlpun_newhis,12,ERR=91) xpos(n),ypos(n),zpos(n)
	    if (keytrj.GT.0) write(dlpun_newhis,12,ERR=91) xvel(n),yvel(n),zvel(n)
	    if (keytrj.GT.1) write(dlpun_newhis,12,ERR=91) xfor(n),yfor(n),zfor(n)
	  end do
	end if
	! Return conditions:
	writeframe=0   ! Success
	goto 99
91	writeframe=-1  ! Failed (error writing file)
99	return
	end function writeframe

	integer function writeheader(fname, un, newfmt)
	implicit none
!	dlpun_newhis = Fortran unit number to read from.
!	newfmt = Format of history file: 0 = Unformatted, 1 = Formatted.
	character*80 :: fname
	integer :: n, un, newfmt
10	FORMAT (A80)
11	FORMAT (2I10)
	his_out_format = newfmt
	dlpun_newhis = un
	if (his_out_format.EQ.-1) stop "HISfile format not yet determined / file not yet opened!"
	if (his_out_format.EQ.0) then   ! Write unformatted history file.....
	  open(unit=dlpun_newhis,file=fname,form="unformatted",status="new")
	  write(dlpun_newhis,ERR=91) dlpheader
	  ! Cell params
	  write(dlpun_newhis,ERR=91) DBLE(natms)
	  ! Now can construct the rest of the atom data part of the HISu file
	  write(dlpun_newhis,ERR=91) (atmname(n),n=1,natms)
	  write(dlpun_newhis,ERR=91) (mass(n),n=1,natms)
	  write(dlpun_newhis,ERR=91) (charge(n),n=1,natms)
	else     ! Write formatted history file....
	  open(unit=dlpun_newhis,file=fname,form="formatted",status="replace")
	  write(dlpun_newhis,10) dlpheader
	  write(dlpun_newhis,11) keytrj,imcon
	end if
	! Return conditions:
	writeheader=0   ! Success
	goto 99
91	writeheader=-1  ! Failed (error writing file)
99	return
	end function writeheader


!	\/\/\/\/\//\/
!	** outinfo **
!	/\/\/\/\/\\/\
!	Subroutine to read molecule type and species info from a DL_POLY OUT file.
	integer function outinfo(fname,verbosity)
	use parse
!	verbosity = Verbosity level: 0 = Quiet, 1 = Full.
	implicit none
	character*80 :: fname,discard,discard2
	integer :: s,discardn,verbosity,n,i

	dlpfn_out = fname
	dlpun_out=19
	maxmols = -1
	maxatoms = -1
	if (dlpfn_out(1:2).EQ."0 ") then  ! Get manual input.
	  write(*,*) "Number of molecule types in config file?"
	  read(*,*) nspecies
          allocate(s_name(nspecies),stat=s); if (s.NE.0) stop "Allocation error <s_name>"
          allocate(s_nmols(nspecies),stat=s); if (s.NE.0) stop "Allocation error <s_nmols>"
          allocate(s_natoms(nspecies),stat=s); if (s.NE.0) stop "Allocation error <s_natoms>"
          allocate(s_start(nspecies),stat=s); if (s.NE.0) stop "Allocation error <s_start>"
          allocate(s_molstart(nspecies),stat=s); if (s.NE.0) stop "Allocation error <s_molstart>"
	  do n=1,nspecies
	    write(*,*) "Number of molecules of type: ",n,"?"
	    read(*,*) s_nmols(n)
	    if (s_nmols(n).GT.maxmols) maxmols = s_nmols(n)
	    write(*,*) "Number of atoms in molecule type: ",n,"?"
	    read(*,*) s_natoms(n)
	    if (s_natoms(n).GT.maxatoms) maxatoms = s_natoms(n)
	  end do
	  goto 75
	end if

	! Read the data from the OUT file....
	open(unit=dlpun_out,file=dlpfn_out,form="formatted",iostat=n)
	if (n.NE.0) then   ! File did not open successfully....
	  write(0,*) "Could not open output file:",dlpfn_out
	  close(dlpun_out)
	  goto 91
	end if
	! Analyse the OUT file: Acquire the molecular data.
	! Allocate arrays as we go along...
50	if (.not.readline(dlpun_out)) goto 91
	if (arg(1).eq."SYSTEM".and.arg(2).eq."SPECIFICATION") goto 51
	!if (arg(1)(1:6).eq."SYSTEM".and.arg(2)(1:13).eq."SPECIFICATION") goto 51
	goto 50
51	if (.not.readline(dlpun_out)) goto 91
	if (arg(1).ne."number".and.arg(3).ne."molecular".and.arg(4).ne."types") goto 51
	nspecies=argi(5)
	if (verbosity.EQ.1) write(0,*) "Number of molecular species in simulation : ",nspecies
	allocate(s_name(nspecies),stat=s); if (s.NE.0) stop "Allocation error <s_name>"
	allocate(s_nmols(nspecies),stat=s); if (s.NE.0) stop "Allocation error <s_nmols>"
	allocate(s_natoms(nspecies),stat=s); if (s.NE.0) stop "Allocation error <s_natoms>"
	allocate(s_start(nspecies),stat=s); if (s.NE.0) stop "Allocation error <s_start>"
	allocate(s_molstart(nspecies),stat=s); if (s.NE.0) stop "Allocation error <s_molstart>"
	do n=1,nspecies
52	  if (.not.readline(dlpun_out)) goto 91
	  if (arg(1).ne."name".and.arg(3).ne."species:") goto 52
	  ! Get the information about this molecular species....
	  s_name(n) = arg(4)
	  do i=5,nargsparsed
	    s_name(n) = trim(s_name(n))//" "//trim(arg(i))
	  end do
	  write(0,*) "Found: ",s_name(n)
53	  if (.not.readline(dlpun_out)) goto 91
	  if (arg(1).ne."number".and.arg(3).ne."molecules") goto 53
	  s_nmols(n)=argi(4)
	  if (verbosity.EQ.1) write(0,*) "-- Number of molecules: ", s_nmols(n)
	  if (s_nmols(n).GT.maxmols) maxmols = s_nmols(n)
54	  if (.not.readline(dlpun_out)) goto 91
	  if (arg(1).ne."number".and.arg(3).ne."atoms/sites") goto 54
	  s_natoms(n)=argi(4)
	  if (verbosity.EQ.1) write(0,*) "-- Number of atoms/sites: ", s_natoms(n)
	end do
	close(dlpun_out)

	! Set the startatoms and startmols for each species.....
75	s_start(1)=1
	s_molstart(1)=1
	do n=2,nspecies
	  s_start(n)=s_start(n-1)+(s_nmols(n-1)*s_natoms(n-1))
	  s_molstart(n)=s_molstart(n-1)+s_nmols(n-1)
	end do

	! Set natoms, maxmols and maxatoms
	natms = 0
	do n=1,nspecies
	  if (s_nmols(n).GT.maxmols) maxmols = s_nmols(n)
	  if (s_natoms(n).GT.maxatoms) maxatoms = s_natoms(n)
	  natms = natms + s_nmols(n)*s_natoms(n)
	end do
	s_totalMols = sum(s_nmols)

	! Allocate other arrays...
	call alloc_xyz
	! Return conditions:
	outinfo=0   ! Success
	goto 99
91	outinfo=-1  ! Failed (error reading file)
99	return
	end function outinfo


!	\/\/\/\/\/\/\/\/\/\/\/
!	** read/writeconfig **
!	/\/\/\/\/\/\/\/\/\/\/\
!	Subroutines to read/write DL_POLY CONFIG files

	subroutine readconfig(fname)
	use parse
	implicit none
	character*80 :: fname
	integer :: n,discardi
	logical :: asframe
	real*8 :: discardr
14      FORMAT (A80)
        ! DL_Poly CONFIG 1 - data in file/image convention/???
16      FORMAT (I10,I10,A100)
        ! DL_POLY CONFIG 2 - atom name/index/atomic number
17      FORMAT (a8,i10,f20.14)
        ! DL_POLY CONFIG 3 - atom positions OR forces OR velocities
19      FORMAT (3f20.14)
	dlpun_cfg = 18
	open(unit=dlpun_cfg,file=fname,form="formatted",status="old",err=100)
        read(dlpun_cfg,14) cfgheader      ! Remove the first line (simulation name)
        read(dlpun_cfg,16) keytrj,imcon,discardi
        ! Check the image convention - if it is NOT zero we read the cell dimensions
        if (imcon.GT.0) then     ! Read in the cell coords if necessary
	  if (.not.readline(dlpun_cfg)) goto 99
	  cfgcell(1) = argr(1); cfgcell(2) = argr(2); cfgcell(3) = argr(3)
	  if (.not.readline(dlpun_cfg)) goto 99
	  cfgcell(4) = argr(1); cfgcell(5) = argr(2); cfgcell(6) = argr(3)
	  if (.not.readline(dlpun_cfg)) goto 99
	  cfgcell(7) = argr(1); cfgcell(8) = argr(2); cfgcell(9) = argr(3)
        end if
        ! Atomic positions/velocities/forces....
	do n=1,natms
	  read(dlpun_cfg,17,err=99) cfgname(n),discardi,discardr
	  read(dlpun_cfg,19,err=99) cfgxyz(n,1:3)
	  if (keytrj.GE.1) read(dlpun_cfg,19,err=99) cfgvel(n,1:3)
	  if (keytrj.EQ.2) read(dlpun_cfg,19,err=99) cfgfrc(n,1:3)
	end do
	close(dlpun_cfg)
	return
99	write(6,"(A,A)") "Error reading file ",fname
	write(6,"(A,I4)") "  Was reading atom ",n
	stop
100	write(6,"(A,A)") "Couldn't open file: ",fname
	stop
	end subroutine readconfig

	subroutine writeconfig(fname)
	implicit none
	character*80 :: fname
	integer :: n
14      FORMAT (A80)
        ! DL_Poly CONFIG 1 - data in file/image convention/???
16      FORMAT (I10,I10,A100)
        ! DL_POLY CONFIG 2 - atom name/index/atomic number
17      FORMAT (A8,2I10)
        ! DL_POLY CONFIG 3 - atom positions OR forces
19      FORMAT (3F20.14)
        ! DL_POLY CONFIG 4 - atom velocities
20      FORMAT (3F20.10)
	dlpun_cfg = 18
	open(unit=dlpun_cfg,file=fname,form="formatted",status="replace",err=100)
        write(dlpun_cfg,14) cfgheader      ! Remove the first line (simulation name)
        write(dlpun_cfg,16) keytrj,imcon
        ! Check the image convention - if it is NOT zero we write the cell dimensions
        if (imcon.gt.0) write(dlpun_cfg,19) cfgcell
	do n=1,natms
	  write(dlpun_cfg,17,err=99) cfgname(n)
	  write(dlpun_cfg,19,err=99) cfgxyz(n,1:3)
	  if (keytrj.GE.1) write(dlpun_cfg,19,err=99) cfgvel(n,1:3)
	  if (keytrj.EQ.2) write(dlpun_cfg,20,err=99) cfgfrc(n,1:3)
	end do
	close(dlpun_cfg)
	return
99	write(6,"(A,A)") "Error while writing file ",fname
	write(6,"(A,I4)") "  Was writing atom ",n
	stop
100	write(6,"(A,A)") "Couldn't write to file: ",fname
	stop
	end subroutine writeconfig

	subroutine alloc_xyz
	implicit none
	integer :: s
	if (xyz_alloc) then
	  write(6,*) "Arrays already allocated!"
	  return
	end if
	allocate(xpos(natms),stat=s); if (s.NE.0) stop "Allocation error for <xpos>"
	allocate(ypos(natms),stat=s); if (s.NE.0) stop "Allocation error for <ypos>"
	allocate(zpos(natms),stat=s); if (s.NE.0) stop "Allocation error for <xpos>"
	allocate(charge(natms),stat=s); if (s.NE.0) stop "Allocation error for <charge>"
	allocate(mass(natms),stat=s); if (s.NE.0) stop "Allocation error for <mass>"
	allocate(xvel(natms),stat=s); if (s.NE.0) stop "Allocation error for <xvel>"
	allocate(yvel(natms),stat=s); if (s.NE.0) stop "Allocation error for <yvel>"
	allocate(zvel(natms),stat=s); if (s.NE.0) stop "Allocation error for <zvel>"
        allocate(xfor(natms),stat=s); if (s.NE.0) stop "Allocation error for <xfor>"
	allocate(yfor(natms),stat=s); if (s.NE.0) stop "Allocation error for <yfor>"
	allocate(zfor(natms),stat=s); if (s.NE.0) stop "Allocation error for <zfor>"
	allocate(atmname(natms),stat=s); if (s.NE.0) stop "Allocation error for <atmnam>"
	! Config arrays
	allocate(cfgxyz(natms,3),stat=s); if (s.NE.0) stop "Allocation error for <cfgxyz>"
	allocate(cfgvel(natms,3),stat=s); if (s.NE.0) stop "Allocation error for <cfgvel>"
	allocate(cfgfrc(natms,3),stat=s); if (s.NE.0) stop "Allocation error for <cfgfrc>"
	allocate(cfgname(natms),stat=s); if (s.NE.0) stop "Allocation error for <cfgname>"
	xyz_alloc = .true.
	end subroutine alloc_xyz

	subroutine dealloc_xyz
	implicit none
	integer :: s
	if (.not.xyz_alloc) then
	  write(6,*) "Arrays not yet allocated."
	  return
	end if
	deallocate(xpos,stat=s); if (s.NE.0) stop "Deallocation error for <xpos>"
	deallocate(ypos,stat=s); if (s.NE.0) stop "Deallocation error for <ypos>"
	deallocate(zpos,stat=s); if (s.NE.0) stop "Deallocation error for <xpos>"
	deallocate(charge,stat=s); if (s.NE.0) stop "Deallocation error for <charge>"
	deallocate(mass,stat=s); if (s.NE.0) stop "Deallocation error for <mass>"
	deallocate(xvel,stat=s); if (s.NE.0) stop "Deallocation error for <xvel>"
	deallocate(yvel,stat=s); if (s.NE.0) stop "Deallocation error for <yvel>"
	deallocate(zvel,stat=s); if (s.NE.0) stop "Deallocation error for <zvel>"
        deallocate(xfor,stat=s); if (s.NE.0) stop "Deallocation error for <xfor>"
	deallocate(yfor,stat=s); if (s.NE.0) stop "Deallocation error for <yfor>"
	deallocate(zfor,stat=s); if (s.NE.0) stop "Deallocation error for <zfor>"
	deallocate(atmname,stat=s); if (s.NE.0) stop "Deallocation error for <atmnam>"
	! Config arrays
	deallocate(cfgxyz,stat=s); if (s.NE.0) stop "Deallocation error for <cfgxyz>"
	deallocate(cfgvel,stat=s); if (s.NE.0) stop "Deallocation error for <cfgvel>"
	deallocate(cfgfrc,stat=s); if (s.NE.0) stop "Deallocation error for <cfgfrc>"
	deallocate(cfgname,stat=s); if (s.NE.0) stop "Deallocation error for <cfgname>"
	xyz_alloc = .false.
	end subroutine dealloc_xyz

!	\/\/\/\/\/\/\/\/\/\/\/
!	** Misc Subroutines **
!	/\/\/\/\/\/\/\/\/\/\/\
!       Other bits and pieces
	character*8 function s_atom(sp,i)
	implicit none
	integer :: sp, i
	s_atom = atmname(s_start(sp)+i-1)
	end function s_atom

	end module dlprw
