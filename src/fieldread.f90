	! Routine to read DL_POLY FIELD file

	module dlpfield

	  integer, parameter :: MAXINTRA = 100, MAXATOMS = 50, MAXRIGID = 3
	  ! Field data (species composition, count etc)
          integer :: f_nspecies, f_totalatoms, f_ntypes
	  integer, allocatable :: f_natoms(:),f_nmols(:),f_start(:),f_nbonds(:),f_nangles(:),f_ntorsions(:),f_nrigid(:)
	  integer, allocatable :: f_bondmatrix(:,:,:), f_nconstraints(:)
	  integer, allocatable :: f_bonds(:,:,:), f_angles(:,:,:), f_torsions(:,:,:), f_types(:,:), f_imasses(:,:)
	  integer, allocatable :: f_rigid(:,:,:), f_hybrid(:,:)
	  real*8, allocatable :: f_masses(:,:)
	  character*8, allocatable :: f_atmname(:,:), f_typenames(:)
	  character*30, allocatable :: f_speciesname(:)
	  real*8, allocatable :: f_charges(:,:), f_vdwscale(:,:,:), f_elecscale(:,:,:)
	  real*8, allocatable :: f_vdwsigma(:,:), f_vdwepsilon(:,:)

	contains

!	\/\/\/\/\/\/\/\/
!	** readfield **
!	/\/\/\/\/\/\/\/\
!	Subroutines to read in (bits of a) DL_POLY FIELD file.

	logical function readfield(fieldfile)
	use parse
	implicit none
	character*8 :: discard
	character*80 :: fieldfile
	integer :: n,m,sp,sp2,i,j, num,last
	logical :: found
	f_totalatoms = 0
	! Open field file
	open(unit=40,file=fieldfile,form="formatted",status="old")
	! Read in and discard title line
	if (.not.readline(40)) goto 999
	! Need to take notice of units line.... otherwise LJ calcs might be wrong (we want it in kJ)
	if (.not.readline(40)) goto 999
	if ((arg(2).ne."kj").AND.(arg(2).ne."kJ")) stop "Energy unit in FIELD file must be in kJ for now"
	! Next line is number of species
	if (.not.readline(40)) goto 999
	if (arg(1).ne."molecules") goto 999
	f_nspecies = argi(2)
	write(6,*) "Number of species defined in FIELD file:", f_nspecies
	allocate(f_natoms(f_nspecies),f_nmols(f_nspecies),f_start(f_nspecies))
	allocate(f_nbonds(f_nspecies),f_nangles(f_nspecies),f_ntorsions(f_nspecies),f_nrigid(f_nspecies))
	allocate(f_types(f_nspecies,MAXATOMS),f_nconstraints(f_nspecies))
	allocate(f_bonds(f_nspecies,MAXINTRA,2),f_angles(f_nspecies,MAXINTRA,3),f_torsions(f_nspecies,MAXINTRA,4))
	allocate(f_rigid(f_nspecies,MAXRIGID,MAXATOMS))
	allocate(f_bondmatrix(f_nspecies,MAXATOMS,MAXATOMS))
	allocate(f_atmname(f_nspecies,MAXATOMS),f_charges(f_nspecies,MAXATOMS),f_speciesname(f_nspecies))
	allocate(f_vdwscale(f_nspecies,MAXATOMS,MAXATOMS),f_elecscale(f_nspecies,MAXATOMS,MAXATOMS))
	allocate(f_typenames(f_nspecies*MAXATOMS), f_imasses(f_nspecies, MAXATOMS), f_masses(f_nspecies, MAXATOMS))
	allocate(f_hybrid(f_nspecies, MAXATOMS))
	f_nbonds = 0
	f_nconstraints = 0
	f_nangles = 0
	f_ntorsions = 0
	f_nrigid = 0
	f_rigid = 0
	f_vdwscale = 1.0
	f_elecscale = 1.0
	f_start = 0
	f_ntypes = 0
	f_types = 0
	f_atmname = "_NOTHING"
	f_masses = 0.0
	f_imasses = 0
	f_bondmatrix = 0
	f_hybrid = 0
	! Loop over species definitions
	do sp=1,f_nspecies
	  ! Name is first
	  if (.not.readline(40)) goto 999
	  f_speciesname(sp) = arg(1)
	  ! Next comes blocks of atoms, bonds, angles, constraints etc.
100	  if (.not.readline(40)) goto 999
	    select case (arg(1))
	      case ("nummols","NUMMOLS")
		f_nmols(sp) = argi(2)
		write(6,*) "...number of molecules of species ",sp,":",f_nmols(sp)
	      case ("atoms","ATOMS")
		f_natoms(sp) = argi(2)
		write(6,*) "...number of atoms in species ",sp,":",f_natoms(sp)
		f_totalatoms = f_totalatoms + f_natoms(sp) * f_nmols(sp)
		if (f_natoms(sp).gt.MAXATOMS) stop "Too many atoms in this species. Increase MAXATOMS."
		do n=1,f_natoms(sp)
		  if (.not.readline(40)) goto 999
		  f_atmname(sp,n) = arg(1)
		  f_masses(sp,n) = argr(2)
		  f_imasses(sp,n) = nint(argr(2))
		  f_charges(sp,n) = argr(3)
		  write(6,*) ".......atomname/imass/charge:", f_atmname(sp,n), f_imasses(sp,n), f_charges(sp,n)
		end do
	      case ("bonds","BONDS")
		write(6,*) "...number of bonds in species ",sp,":",argi(2)
		do n=1,argi(2)
		  f_nbonds(sp) = f_nbonds(sp) + 1
		  if (f_nbonds(sp).gt.MAXINTRA) stop "Too many bonds in this species. Increase MAXINTRA."
		  if (.not.readline(40)) goto 999
		  f_bonds(sp,f_nbonds(sp),1) = argi(2)
		  f_bonds(sp,f_nbonds(sp),2) = argi(3)
		  f_bondmatrix(sp,f_bonds(sp,f_nbonds(sp),1),f_bonds(sp,f_nbonds(sp),2)) = 1
		  f_bondmatrix(sp,f_bonds(sp,f_nbonds(sp),2),f_bonds(sp,f_nbonds(sp),1)) = 1
		  write(6,"(a,4i5)") ".......bond i-j:", f_bonds(sp,f_nbonds(sp),1), f_bonds(sp,f_nbonds(sp),2)
		end do
	      case ("constraints","CONSTRAINTS")
		write(6,*) "...number of constraints in species ",sp,":",argi(2)
		do n=1,argi(2)
		  f_nbonds(sp) = f_nbonds(sp) + 1
		  f_nconstraints(sp) = f_nconstraints(sp) + 1
		  if (f_nbonds(sp).gt.MAXINTRA) stop "Too many bonds in this species. Increase MAXINTRA."
		  if (.not.readline(40)) goto 999
		  f_bonds(sp,f_nbonds(sp),1) = argi(1)
		  f_bonds(sp,f_nbonds(sp),2) = argi(2)
		  f_bondmatrix(sp,f_bonds(sp,f_nbonds(sp),1),f_bonds(sp,f_nbonds(sp),2)) = 1
		  f_bondmatrix(sp,f_bonds(sp,f_nbonds(sp),2),f_bonds(sp,f_nbonds(sp),1)) = 1
		  write(6,"(a,4i5)") ".......constraint i-j:", f_bonds(sp,f_nbonds(sp),1), f_bonds(sp,f_nbonds(sp),2)
		end do
	      case ("angles","ANGLES")
		write(6,*) "...number of angles in species ",sp,":",argi(2)
		do n=1,argi(2)
		  f_nangles(sp) = f_nangles(sp) + 1
		  if (f_nangles(sp).gt.MAXINTRA) stop "Too many angles in this species. Increase MAXINTRA."
		  if (.not.readline(40)) goto 999
		  f_angles(sp,f_nangles(sp),1) = argi(2)
		  f_angles(sp,f_nangles(sp),2) = argi(3)
		  f_angles(sp,f_nangles(sp),3) = argi(4)
		  write(6,"(a,4i5)") ".......angle i-j-k:", (f_angles(sp,f_nangles(sp),m), m=1,3)
		end do
	      case ("dihedrals","DIHEDRALS")
		write(6,*) "...number of torsions in species ",sp,":",argi(2)
		do n=1,argi(2)
		  f_ntorsions(sp) = f_ntorsions(sp) + 1
		  if (f_ntorsions(sp).gt.MAXINTRA) stop "Too many torsions in this species. Increase MAXINTRA."
		  if (.not.readline(40)) goto 999
		  f_torsions(sp,f_ntorsions(sp),1) = argi(2)
		  f_torsions(sp,f_ntorsions(sp),2) = argi(3)
		  f_torsions(sp,f_ntorsions(sp),3) = argi(4)
		  f_torsions(sp,f_ntorsions(sp),4) = argi(5)
		  write(6,"(a,4i5)") ".......torsion i-j-k-l:", (f_torsions(sp,f_ntorsions(sp),m), m=1,4)
		  ! Apply the scaling factors to the matrices now...
		  f_elecscale(sp,f_torsions(sp,f_ntorsions(sp),1),f_torsions(sp,f_ntorsions(sp),4)) = argr(9)
		  f_vdwscale(sp,f_torsions(sp,f_ntorsions(sp),1),f_torsions(sp,f_ntorsions(sp),4)) = argr(10)
		  f_elecscale(sp,f_torsions(sp,f_ntorsions(sp),4),f_torsions(sp,f_ntorsions(sp),1)) = argr(9)
		  f_vdwscale(sp,f_torsions(sp,f_ntorsions(sp),4),f_torsions(sp,f_ntorsions(sp),1)) = argr(10)
		end do
	      case ("rigid","RIGID")
		do n=1,argi(2)
		  f_nrigid(sp) = f_nrigid(sp) + 1
		  if (.not.readline(40)) goto 999
		  ! First number is the number of atoms in this rigid unit
		  i = 2
		  num = argi(1)
		  do m=1,num
		     f_rigid(sp,f_nrigid(sp),m) = argi(i)
		     i = i + 1
		     ! Read in a new line if we've reached the end of this one (16th item)
		     if ((mod(m,16).eq.0).and.(m.lt.num)) then
			i = 1
			if (.not.readline(40)) goto 999
		     end if
		  end do
		end do
	      case ("finish","FINISH")
		if (sp.gt.1) f_start(sp) = f_start(sp-1) + f_natoms(sp-1)*f_nmols(sp-1)
		cycle
	      case default
		write(6,*) "...unrecognised FIELD entry ", arg(1)
		stop
	    end select
	  goto 100
	end do
	! VDW keyword should follow, so find and assign unique type names/ids
	f_ntypes = 0
	
	do sp=1,f_nspecies
	  do i=1,f_natoms(sp)
	    !write(6,*) sp,i
	    ! Check to see if this atom name occurs before this occurrence
	    found = .FALSE.
	    do sp2=1,sp
	      last = f_natoms(sp2)
	      if (sp.eq.sp2) last = i-1
	      do j=1,last
	        if (f_atmname(sp,i).eq.f_atmname(sp2,j)) found = .TRUE.
	        if (found) exit
	      end do
	      if (found) exit
	    end do
	    ! Is this the first occurrence?
	    if (found) then
	      f_types(sp,i) = f_types(sp2,j)
	    else
	      f_ntypes = f_ntypes + 1
	      f_typenames(f_ntypes) = f_atmname(sp,i)
	      f_types(sp,i) = f_ntypes
	    end if
	  end do
	end do

! 	do n=1,MAXATOMS*f_nspecies
! 	  sp=n/(MAXATOMS+1)+1
! 	  i=mod(n-1,MAXATOMS)+1
! 	  write(6,*) sp,i
! 	  if (f_atmname(sp,i).eq."_NOTHING") cycle
! 	  ! Check to see if this atom name occurs before this occurrence
! 	  found = .FALSE.
! 	  do m=1,n-1
! 	    sp2=m/(MAXATOMS+1)+1
! 	    j=mod(m-1,MAXATOMS)+1
! 	    if (f_atmname(sp2,j).eq."_NOTHING") cycle
! 	    if (f_atmname(sp,i).eq.f_atmname(sp2,j)) found = .TRUE.
! 	    if (found) exit
! 	  end do
! 	  ! Is this the first occurrence?
! 	  if (found) then
! 	    f_types(sp,i) = f_types(sp2,j)
! 	  else
! 	    f_ntypes = f_ntypes + 1
! 	    f_typenames(f_ntypes) = f_atmname(sp,i)
! 	    f_types(sp,i) = f_ntypes
! 	  end if
! 	end do
	write(6,*) "Number of unique atom type (names) found:",f_ntypes
! 	do sp=1,f_nspecies
! 	  do i=1,f_natoms(sp)
! 	   write(6,*) sp,i,f_types(sp,i)
! 	  end do
! 	end do
	allocate(f_vdwsigma(f_ntypes,f_ntypes),f_vdwepsilon(f_ntypes,f_ntypes))
	f_vdwsigma = 0.0d0
	f_vdwepsilon = 0.0d0
	if (.not.readline(40)) goto 999
	if ((arg(1).ne."vdw").and.(arg(1).ne."VDW").and.(arg(1).ne."nvdw").and.(arg(1).ne."NVDW")) goto 999
	do n=1,argi(2)
	  if (.not.readline(40)) goto 999
	  i = 0
	  j = 0
	  do m=1,f_ntypes
	    if (f_typenames(m).eq.arg(1)) then
	      i = m
	      exit
	    end if
	  end do
	  do m=1,f_ntypes
	    if (f_typenames(m).eq.arg(2)) then
	      j = m
	      exit
	    end if
	  end do
	  ! Sanity check - did we find both types in our own list?
	  if ((i.eq.0).or.(j.eq.0)) stop "Error - failed to find type in type list when reading VDW"
	  f_vdwepsilon(i,j) = argr(4)
	  f_vdwepsilon(j,i) = argr(4)
	  f_vdwsigma(i,j) = argr(5)
	  f_vdwsigma(j,i) = argr(5)
	end do

	close(40)

	! Finalise scaling matrices - torsions have already been processed, every other array element will be 1.0
	do sp=1,f_nspecies
	  do n=1,f_nbonds(sp)
	    f_elecscale(sp,f_bonds(sp,n,1),f_bonds(sp,n,2)) = 0.0
	    f_vdwscale(sp,f_bonds(sp,n,1),f_bonds(sp,n,2)) = 0.0
	    f_elecscale(sp,f_bonds(sp,n,2),f_bonds(sp,n,1)) = 0.0
	    f_vdwscale(sp,f_bonds(sp,n,2),f_bonds(sp,n,1)) = 0.0
	  end do
	  do n=1,f_nangles(sp)
	    f_elecscale(sp,f_angles(sp,n,1),f_angles(sp,n,3)) = 0.0
	    f_vdwscale(sp,f_angles(sp,n,1),f_angles(sp,n,3)) = 0.0
	    f_elecscale(sp,f_angles(sp,n,3),f_angles(sp,n,1)) = 0.0
	    f_vdwscale(sp,f_angles(sp,n,3),f_angles(sp,n,1)) = 0.0
	  end do
	  do n=1,f_nrigid(sp)
	    do i=1,MAXATOMS
	      if (f_rigid(sp,n,i).eq.0) exit
	      do j=i+1,MAXATOMS
		if (f_rigid(sp,n,j).eq.0) exit
		f_elecscale(sp,i,j) = 0.0
		f_vdwscale(sp,i,j) = 0.0
		f_elecscale(sp,j,i) = 0.0
		f_vdwscale(sp,j,i) = 0.0
	      end do
	    end do
	  end do
	end do
	
! 	do n=1,6
! 	  write(6,"(6f8.4)") (f_elecscale(1,n,m),m=1,6)
! 	end do

	readfield = .TRUE.
	goto 1000
999	write(6,*) "Error reading FIELD file"
	readfield = .FALSE.
1000	return
	end function readfield

	integer function getpartner(sp, i, notj)
	implicit none
	integer :: n, sp, i, notj
	getpartner = 0
	do n=1,f_natoms(sp)
	  if ((f_bondmatrix(sp,i,n).eq.1).and.(n.ne.notj)) then
	    getpartner = n
	    return
	  end if
	end do
	end function getpartner

	end module dlpfield
