	! Useful routines for many other programs!

	module utility
	  ! Centre-of-mass arrays
	  real*8, allocatable :: comx(:,:),comy(:,:),comz(:,:)
	  ! Axis Arrays
	  integer, allocatable :: axesAatoms(:,:), axesBatoms(:,:)
	  logical, allocatable :: axesAdefined(:), axesBdefined(:)
	  real*8, allocatable :: axesA(:,:,:), axesB(:,:,:)
	  real*8, allocatable :: axesAorigin(:,:,:), axesBorigin(:,:,:)
	  logical :: com_alloc = .false., axis_alloc = .false.
	  ! Reciprocal cell
	  real*8 :: rcell(9), celldet
	  ! Inverse cell
	  real*8 :: icell(9)
	contains

	! ========================================
	! Text / String Handling
	! ========================================
	! Return constructed output filename
	character*80 function outputFileName(baseName, defaultName, suffix)
	implicit none
	character*80, intent(in) :: baseName
	character*(*), intent(in) :: defaultName, suffix
	integer :: baseLength, n, dotPos

	if (LEN_TRIM(baseName).eq.0) then
	  outputFileName = defaultName
	else if (baseName.eq."0") then
	  outputFileName = defaultName
	else
	  dotPos = SCAN(baseName,".",.true.)
	  if (dotPos.eq.0) then
	    outputFileName = TRIM(baseName)
	  else
	    outputFileName = baseName(1:dotPos-1)
	  end if
	end if

	! Add on trailing dot
	outputFileName = TRIM(outputFileName)//"."
	baseLength = LEN_TRIM(outputFileName)

	outputFileName = TRIM(outputFileName)//TRIM(suffix)
	end function outputFileName

	! Return formatted string of the form N_MM
	character*4 function stringNMM(n, m)
	implicit none
	integer, intent(in) :: n, m

	if ((n.gt.9).or.(n.lt.0)) write(0,*) "WARNING : stringNMM - 'n' is out of range (0-9)"
	if ((m.gt.99).or.(m.lt.0)) write(0,*) "WARNING : stringNMM - 'm' is out of range (0-99)"
	stringNMM(1:1) = CHAR(48+n)
	stringNMM(2:2) = '_'
	stringNMM(3:3) = CHAR(48+(m/10))
	stringNMM(4:4) = CHAR(48+mod(m,10))
	end function stringNMM

	! Return formatted string of the form N_MM_OO
	character*7 function stringNMMOO(n, m, o)
	implicit none
	integer, intent(in) :: n, m, o

	if ((n.gt.9).or.(n.lt.0)) write(0,*) "WARNING : stringNMMOO - 'n' is out of range (0-9)"
	if ((m.gt.99).or.(m.lt.0)) write(0,*) "WARNING : stringNMMOO - 'm' is out of range (0-99)"
	if ((o.gt.99).or.(o.lt.0)) write(0,*) "WARNING : stringNMMOO - 'o' is out of range (0-99)"
	stringNMMOO(1:1) = CHAR(48+n)
	stringNMMOO(2:2) = '_'
	stringNMMOO(3:3) = CHAR(48+(m/10))
	stringNMMOO(4:4) = CHAR(48+mod(m,10))
	stringNMMOO(5:5) = '_'
	stringNMMOO(6:6) = CHAR(48+(o/10))
	stringNMMOO(7:7) = CHAR(48+mod(o,10))
	end function stringNMMOO

	! ========================================
	! Cell / Reciprocal Cell
	! ========================================
	subroutine calc_rcell()
	use dlprw; implicit none
	real*8 :: rvol
	real*8, parameter :: pi = 3.14159265358979d0
	integer :: n
	if (imcon.gt.3) stop "calc_rcell can't handle this image convention."

	! Calculate the reciprocal cell (cubic, orthorhombic, or parallelepiped)
	! Take cross products of original cell vectors to get reciprocal cell vectors
	! R(1,4,7) = C(Y) * C(Z)
	! R(2,5,8) = C(Z) * C(X)
	! R(3,6,9) = C(X) * C(Y)
	rcell(1) = cell(5) * cell(9) - cell(8) * cell(6)
	rcell(4) = cell(8) * cell(3) - cell(2) * cell(9)
	rcell(7) = cell(2) * cell(6) - cell(5) * cell(3)

	rcell(2) = cell(6) * cell(7) - cell(9) * cell(4)
	rcell(5) = cell(9) * cell(1) - cell(3) * cell(7)
	rcell(8) = cell(3) * cell(4) - cell(6) * cell(1)

	rcell(3) = cell(4) * cell(8) - cell(7) * cell(5)
	rcell(6) = cell(7) * cell(2) - cell(1) * cell(8)
	rcell(9) = cell(1) * cell(5) - cell(4) * cell(2)
	rvol = 2.0d0*pi / volume(cell)
	rcell = rcell * rvol
	write(0,"(3F12.4)") rcell
	end subroutine calc_rcell

	subroutine calc_icell()
	use dlprw; implicit none
	if (has_icell) return
	icell = cell
	call gaussj(icell,3,3)
	has_icell = .true.
	end subroutine calc_icell

	real*8 function volume(cell)
	implicit none
	real*8, intent(in) :: cell(9)
	! Hard-coded calculation of determinant of 3x3 matrix
	volume = cell(1) * (cell(5)*cell(9) - cell(8)*cell(6));
	volume = volume - cell(2) * (cell(4)*cell(9) - cell(7)*cell(6));
	volume = volume + cell(3) * (cell(4)*cell(8) - cell(7)*cell(5));
	end function volume

	! ========================================
	! Periodic boundary conditions
	! ========================================
	subroutine pbc(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	use dlprw
	implicit none
	! Performs minimum image convention.
	! xyz1: point to consider; xyz2: reference point; xyz3: minimum image coords of xyz1 (result)
	real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x,y,z
	! Performs minimum image convention.....
	if (imcon.eq.0) then
	  x3=x1
	  y3=y1
	  z3=z1
	else if ((imcon.eq.1).or.(imcon.eq.2)) then
	  x3=x1 - cell(1)*NINT((x1-x2)/cell(1))
	  y3=y1 - cell(5)*NINT((y1-y2)/cell(5))
	  z3=z1 - cell(9)*NINT((z1-z2)/cell(9))
	else if (imcon.eq.3) then
	  ! Make sure we have an up-to-date inverse cell
	  call calc_icell()
	  ! Multiply coordinate 
	  !write(0,"(a,3f12.4)") "Original coords = ",x1,y1,z1
	  x = x1 - x2
	  y = y1 - y2
	  z = z1 - z2
	  call vmatmult(x,y,z,icell)
	  if (x.lt.-0.5) x = x + 1.0
	  if (x.gt.0.5) x = x - 1.0
	  if (y.lt.-0.5) y = y + 1.0
	  if (y.gt.0.5) y = y - 1.0
	  if (z.lt.-0.5) z = z + 1.0
	  if (z.gt.0.5) z = z - 1.0
	  call vmatmult(x,y,z,cell)
	  x3 = x + x2
	  y3 = y + y2
	  z3 = z + z2
	  !write(0,"(a,3f12.4)") "MIM coords = ",x3,y3,z3
	else
	  write(0,*) "This boundary condition is not supported."
	end if
	end subroutine pbc

	subroutine dlpfold(x1,y1,z1,x2,y2,z2)
	use dlprw
	implicit none
	! Fold coordinates into DL_POLY unit cell (-0.5 -> 0.5)
	! xyz1: point to fold; xyz2: folded coordinate (result)
	real*8 :: x1,y1,z1,x2,y2,z2
	if (imcon.eq.0) then
	  x2=x1
	  y2=y1
	  z2=z1
	else if ((imcon.eq.1).or.(imcon.eq.2)) then
	  x2=x1
	  y2=y1
	  z2=z1
	  if (x2.lt.-0.5*cell(1)) x2 = x2 + cell(1)
	  if (x2.gt.0.5*cell(1)) x2 = x2 - cell(1)
	  if (y2.lt.-0.5*cell(5)) y2 = y2 + cell(5)
	  if (y2.gt.0.5*cell(5)) y2 = y2 - cell(5)
	  if (z2.lt.-0.5*cell(9)) z2 = z2 + cell(9)
	  if (z2.gt.0.5*cell(9)) z2 = z2 - cell(9)
	else if (imcon.eq.3) then
	  ! Make sure we have an up-to-date inverse cell
	  call calc_icell()
	  ! Multiply coordinate 
	  !write(0,"(a,3f12.4)") "Original coords = ",x1,y1,z1
	  x2 = x1
	  y2 = y1
	  z2 = z1
	  call vmatmult(x2,y2,z2,icell)
	  if (x2.lt.-0.5) x2 = x2 + 1.0
	  if (x2.gt.0.5) x2 = x2 - 1.0
	  if (y2.lt.-0.5) y2 = y2 + 1.0
	  if (y2.gt.0.5) y2 = y2 - 1.0
	  if (z2.lt.-0.5) z2 = z2 + 1.0
	  if (z2.gt.0.5) z2 = z2 - 1.0
	  call vmatmult(x2,y2,z2,cell)
	else
	  write(0,*) "This boundary condition is not supported."
	end if
	end subroutine dlpfold

	subroutine unitfold(x1,y1,z1)
	use dlprw
	implicit none
	! Fold coordinates into unit cell (0 -> 1)
	! xyz1: point to fold
	real*8 :: x1,y1,z1
	if (imcon.eq.0) then
	  ! No cell, do nothing
	else if ((imcon.eq.1).or.(imcon.eq.2)) then
	  if (x1.lt.0.0) x1 = x1 - int(x1/cell(1)-1)*cell(1)
	  if (x1.gt.cell(1)) x1 = x1 - int(x1/cell(1)+1)*cell(1)
	  if (y1.lt.0.0) y1 = y1 - int(y1/cell(5)-1)*cell(5)
	  if (y1.gt.cell(5)) y1 = y1 - int(y1/cell(5)+1)*cell(5)
	  if (z1.lt.0.0) z1 = z1 - int(z1/cell(9)-1)*cell(9)
	  if (z1.gt.cell(9)) z1 = z1 - int(z1/cell(9)+1)*cell(9)
	!else if (imcon.eq.3) then
	else
	  write(0,*) "This boundary condition is not supported."
	end if
	end subroutine unitfold

	! ========================================
	! Centre of mass (arrays)
	! ========================================
	subroutine calc_com
	use dlprw
	implicit none
	! Calculate the centre of mass for each molecule of all species
	integer :: sp
	if (.NOT.com_alloc) call alloc_com(nspecies,maxval(s_nmols))
	do sp=1,nspecies
	  call calc_spcom(sp)
	end do
	end subroutine calc_com

	subroutine calc_spcom(sp)
	use dlprw
	implicit none
	! Calculate the centre of mass for each molecule of all species
	real*8 :: massnorm,tx,ty,tz,x0,y0,z0,tempmass
	integer :: sp,n,m,i
	if (.NOT.com_alloc) call alloc_com(nspecies,maxval(s_nmols))
	i = s_start(sp)
	do m=1,s_nmols(sp)
	  ! Calculate COM relative to first atom in molecule
	  tempmass = mass(i)
	  if (mass(i).lt.1.0e-3) tempmass = 1.0
	  comx(sp,m) = xpos(i)*tempmass
	  comy(sp,m) = ypos(i)*tempmass
	  comz(sp,m) = zpos(i)*tempmass
	  massnorm = tempmass
	  if (mass(i).lt.1.0e-3) massnorm = 1.0
	!write(66,"(2i6,3f12.6)") m,1,comx(sp,m),comy(sp,m),comz(sp,m)
	  do n=i+1,i+s_natoms(sp)-1
	    call PBC(xpos(n),ypos(n),zpos(n),xpos(i),ypos(i),zpos(i),tx,ty,tz)
	!write(66,"(2i6,3f12.6)") m,n,tx,ty,tz
	    tempmass = mass(n)
	    if (mass(n).lt.1.0e-3) tempmass = 1.0
	    comx(sp,m) = comx(sp,m) + tx*tempmass
	    comy(sp,m) = comy(sp,m) + ty*tempmass
	    comz(sp,m) = comz(sp,m) + tz*tempmass
	    massnorm = massnorm + tempmass
	  end do
	  comx(sp,m) = comx(sp,m) / massnorm
	  comy(sp,m) = comy(sp,m) / massnorm
	  comz(sp,m) = comz(sp,m) / massnorm
	!write(66,"(a,3f12.6)") "com",comx(sp,m),comy(sp,m),comz(sp,m) 
	!write(66,"(2i6,3f12.6)") m,1,comx(sp,m),comy(sp,m),comz(sp,m)
	  i = i + s_natoms(sp)
	end do
	end subroutine calc_spcom

	subroutine alloc_com(i,j)
	implicit none; integer :: n,i,j,status
	com_alloc = .true.
	allocate(comx(i,j),stat=status); if (status.GT.0) stop "Allocation error for comx()"
	allocate(comy(i,j),stat=status); if (status.GT.0) stop "Allocation error for comy()"
	allocate(comz(i,j),stat=status); if (status.GT.0) stop "Allocation error for comz()"
	end subroutine alloc_com

	! ========================================
	! Centre of geometry
	! ========================================
	! Calculate average, mim'd coordinates from array of atom indices provided
	subroutine averagePosition(array, nItems, offset, v)
	use dlprw
	implicit none
	integer, intent(in) :: nItems, array(nItems)
	integer, intent(in) :: offset
	integer :: n, i
	real*8, intent(out) :: v(3)
	real*8 :: t(3), ref(3)

	v = 0.0
	if (nItems.eq.0) return

	! Calculate COG relative to first atom in list
	ref(1) = xpos(array(1)+offset)
	ref(2) = ypos(array(1)+offset)
	ref(3) = zpos(array(1)+offset)
	v(1) = ref(1)
	v(2) = ref(2)
	v(3) = ref(3)
	do n=2,nItems
	  i = offset + array(n)
	  call PBC(xpos(i),ypos(i),zpos(i),ref(1),ref(2),ref(3),t(1),t(2),t(3))
	  v = v + t
	end do

	v = v / real(nItems)

	end subroutine averagePosition

	! ========================================
	! Calculate Axes
	! ========================================
	subroutine calculate_axes(aa, axes, origin)
	use dlprw
	implicit none
	integer, intent(in) :: aa(4)
	real*8, intent(out) :: axes(9), origin(3)
	real*8 :: dp, tx, ty, tz, vx, vy, vz, xmag, xmagsq, ymag
	! 1) Determine the X axis of the molecule
	!    -- Get the vector between atoms 1 and 2...
	call pbc(xpos(aa(2)), ypos(aa(2)), zpos(aa(2)), xpos(aa(1)), ypos(aa(1)), zpos(aa(1)), tx, ty, tz)
	call getVector(tx,ty,tz, xpos(aa(1)), ypos(aa(1)), zpos(aa(1)), vx, vy, vz)
	xmag = SQRT(vx**2 + vy**2 + vz**2)     ! Magnitude of vector
	xmagsq = xmag * xmag
	axes(1) = vx		! Store the un-normalised components of the x-axis vector
	axes(2) = vy
	axes(3) = vz
	! 2) Set the origin of the axis system - halfway along the x-axis vector
	origin(1) = xpos(aa(1)) + 0.5*vx		! N.B. *NOT* the normalised axis, but the actual distance i->j
	origin(2) = ypos(aa(1)) + 0.5*vy
	origin(3) = zpos(aa(1)) + 0.5*vz
	! 3) Determine the Y axis of the molecule
	!    -- Get the vector between midway along the x-axis vector and midway between
	!	the positions of aa 3 and 4.
	call pbc(xpos(aa(3)), ypos(aa(3)), zpos(aa(3)), origin(1), origin(2), origin(3), tx, ty, tz)
	call pbc(xpos(aa(4)), ypos(aa(4)), zpos(aa(4)), origin(1), origin(2), origin(3), vx, vy, vz)	! Use v_xyz temporarily
	tx = 0.5*(tx+vx)
	ty = 0.5*(ty+vy)
	tz = 0.5*(tz+vz)
	call getVector(origin(1),origin(2),origin(3),tx,ty,tz,vx,vy,vz)
	axes(4) = vx   ! Store the un-normalised components of the y-axis vector
	axes(5) = vy
	axes(6) = vz
	! 3a) Orthogonalise the y vector w.r.t. the x vector via Gram-Schmidt
	dp = axes(1)*axes(4) + axes(2)*axes(5) + axes(3)*axes(6)
	tx = axes(4) - (dp / xmagsq)*axes(1)
	ty = axes(5) - (dp / xmagsq)*axes(2)
	tz = axes(6) - (dp / xmagsq)*axes(3)
	axes(4) = tx
	axes(5) = ty
	axes(6) = tz
	ymag=SQRT(tx**2 + ty**2 + tz**2)     ! Magnitude of vector
	! 3b) Normalise the y and x vectors
	axes(1:3) = axes(1:3) / xmag
	axes(4:6) = axes(4:6) / ymag
	! 4) Calculate the z axis from the x and y axes....
	axes(7) = axes(2)*axes(6) - axes(3)*axes(5)
	axes(8) = axes(3)*axes(4) - axes(1)*axes(6)
	axes(9) = axes(1)*axes(5) - axes(2)*axes(4)
	!write(88,*) "O: ",axisox(sp,n),axisoy(sp,n),axisoz(sp,n)
	!write(88,*) "X: ",axisx(sp,n,1),axisx(sp,n,2),axisx(sp,n,3)
	!write(88,*) "Y: ",axisy(sp,n,1),axisy(sp,n,2),axisy(sp,n,3)
	!write(88,*) "Z: ",axisz(sp,n,1),axisz(sp,n,2),axisz(sp,n,3)
	end subroutine calculate_axes

	subroutine genaxes
	use dlprw
	implicit none
	! Subroutine to calculate individual molecular axes of species molecules
	real*8 :: vx,vy,vz,tx,ty,tz
	real*8 :: xmag, ymag, dp
	integer :: sp,m,i
	if (.NOT.axis_alloc) stop "Axes not allocated, but genaxes() was called."
	do sp=1,nspecies
	  ! Calculate local axes if we can
	  ! -- Axes A
	  if (axesAdefined(sp).and.(s_natoms(sp).gt.2).and.(axesAatoms(sp,1).ne.-1)) then
	    i = s_start(sp) - 1
	    do m=1,s_nmols(sp)
	      call calculate_axes(axesAatoms(sp,:)+i, axesA(sp,m,:), axesAorigin(sp,m,:))
	      i = i + s_natoms(sp)
	    end do
	  else
	    ! No axis definition (or cartesian requested) - Set the axes to be in standard orientation
	    call calc_spcom(sp)
	    do m=1,s_nmols(sp)
	      axesA(sp,m,1:3) = (/ 1.0,0.0,0.0 /)
	      axesA(sp,m,4:6) = (/ 0.0,1.0,0.0 /)
	      axesA(sp,m,7:9) = (/ 0.0,0.0,1.0 /)
	      axesAorigin(sp,m,1) = comx(sp,m)
	      axesAorigin(sp,m,2) = comy(sp,m)
	      axesAorigin(sp,m,3) = comz(sp,m)
	    end do
	  end if
	  ! -- Axes B
	  if (axesBdefined(sp).and.(s_natoms(sp).gt.2)) then
	    i = s_start(sp) - 1
	    do m=1,s_nmols(sp)
	      call calculate_axes(axesBatoms(sp,:)+i, axesB(sp,m,:), axesBorigin(sp,m,:))
	      i = i + s_natoms(sp)
	    end do
	  end if
	end do
	end subroutine genaxes

	subroutine alloc_axis()
	use dlprw; implicit none
	integer :: status
	axis_alloc = .true.
	allocate(axesA(nspecies,maxmols,9),stat=status); if (status.GT.0) stop "Allocation error for axisA()"
	allocate(axesAorigin(nspecies,maxmols,3),stat=status); if (status.GT.0) stop "Allocation error for axisAorigin()"
	allocate(axesAatoms(nspecies,4))
	allocate(axesAdefined(nspecies))
	allocate(axesB(nspecies,maxmols,9),stat=status); if (status.GT.0) stop "Allocation error for axisB()"
	allocate(axesBorigin(nspecies,maxmols,3),stat=status); if (status.GT.0) stop "Allocation error for axisBorigin()"
	allocate(axesBatoms(nspecies,4))
	allocate(axesBdefined(nspecies))
	axesAdefined = .false.
	axesBdefined = .false.
	axesAatoms = 0
	axesBatoms = 0
	end subroutine alloc_axis

	! ========================================
	! 3-Vectors
	! ========================================

	subroutine getVector(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	implicit none
	real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3
	! xyz1,xyz2: coordinates of two atoms; xyz3: vector between them (result)
	x3=x1-x2
	y3=y1-y2
	z3=z1-z2
	end subroutine getVector

	real*8 function vec3Magnitude(vec)
	implicit none
	real*8 :: vec(3)
	vec3Magnitude = dsqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
	end function vec3Magnitude
     
	real*8 function vec3DotProduct(abc,xyz)
	implicit none
	real*8 :: abc(3), xyz(3), result
	result=abc(1)*xyz(1) + abc(2)*xyz(2) + abc(3)*xyz(3)
	vec3DotProduct = result
	end function vec3DotProduct

	subroutine vec3CrossProduct(abc,xyz,result)
	implicit none
	real*8 :: abc(3), xyz(3), result(3)
	result(1)=abc(2)*xyz(3) - abc(3)*xyz(2)
	result(2)=abc(3)*xyz(1) - abc(1)*xyz(3)
	result(3)=abc(1)*xyz(2) - abc(2)*xyz(1)
	end subroutine vec3CrossProduct

	subroutine vmatmult(x1,y1,z1,a)
	implicit none
	real*8 :: x1,y1,z1,x2,y2,z2,a(9)
	! Multiply vector(3) v by matrix(3x3) a
	! Assume matrix is stored in row-major order
	x2 = x1*a(1) + y1*a(4) + z1*a(7)
	y2 = x1*a(2) + y1*a(5) + z1*a(8)
	z2 = x1*a(3) + y1*a(6) + z1*a(9)
	x1 = x2
	y1 = y2
	z1 = z2
	end subroutine vmatmult

	! ========================================
	! Geometry
	! ========================================
	real*8 function calculateTorsion(ii,jj,kk,ll,offset)
	use dlprw
	implicit none
	real*8, parameter :: radcon = 57.29577951d0
	integer, intent(in) :: ii, jj, kk, ll, offset
	integer :: i, j, k, l
	real*8 :: tx, ty, tz, ktx, kty, ktz, dp, mag1, mag2
	real*8 :: vecji(3), vecjk(3), veckj(3), veckl(3), xp1(3), xp2(3)

	! Minimum Image w.r.t. atom j 
	i = ii + offset
	j = jj + offset
	k = kk + offset
	l = ll + offset

	! Angle i-j-k
	! Vector j->i
	call pbc(xpos(i),ypos(i),zpos(i),xpos(j),ypos(j),zpos(j),tx,ty,tz)
	call getVector(tx,ty,tz,xpos(j),ypos(j),zpos(j),vecji(1),vecji(2),vecji(3))
	! Vector j->k
	call pbc(xpos(k),ypos(k),zpos(k),xpos(j),ypos(j),zpos(j),ktx,kty,ktz)
	call getVector(ktx,kty,ktz,xpos(j),ypos(j),zpos(j),vecjk(1),vecjk(2),vecjk(3))
	! Angle j-k-l (mim w.r.t. k (mim j))
	! Vector k->j
	veckj = -vecjk
	! Vector k->l
	call pbc(xpos(l),ypos(l),zpos(l),ktx,kty,ktz,tx,ty,tz)
	call getVector(tx,ty,tz,ktx,kty,ktz,veckl(1),veckl(2),veckl(3))

	! Calculate cross products and magnitudes
	call vec3CrossProduct(vecji,vecjk,xp1)
	mag1 = vec3Magnitude(xp1)
	call vec3CrossProduct(veckj,veckl,xp2)
	mag2 = vec3Magnitude(xp2)
		  
	! Calculate dot product and angle...
	dp = vec3DotProduct(xp1, xp2) / (mag1 * mag2)
	calculateTorsion = acos(dp)*radcon

	! Calculate sign
	dp = vec3DotProduct(xp1, veckl)
	if (dp.lt.0) calculateTorsion = -calculateTorsion

	end function calculateTorsion

	real*8 function calculateAngle(ii,jj,kk,offset)
	use dlprw
	implicit none
	real*8, parameter :: radcon = 57.29577951d0
	integer, intent(in) :: ii, jj, kk, offset
	integer :: i, j, k
	real*8 :: tx, ty, tz, ktx, kty, ktz, dp, mag1, mag2
	real*8 :: vecji(3), vecjk(3)

	! Minimum Image w.r.t. atom j 
	i = ii + offset
	j = jj + offset
	k = kk + offset

	! Vector j->i
	call pbc(xpos(i),ypos(i),zpos(i),xpos(j),ypos(j),zpos(j),tx,ty,tz)
	call getVector(tx,ty,tz,xpos(j),ypos(j),zpos(j),vecji(1),vecji(2),vecji(3))
	mag1 = vec3Magnitude(vecji)
	! Vector j->k
	call pbc(xpos(k),ypos(k),zpos(k),xpos(j),ypos(j),zpos(j),tx,ty,tz)
	call getVector(tx,ty,tz,xpos(j),ypos(j),zpos(j),vecjk(1),vecjk(2),vecjk(3))
	mag2 = vec3Magnitude(vecjk)

	! Calculate dot product and angle...
	dp = vec3DotProduct(vecji, vecjk) / (mag1 * mag2)
	calculateAngle = acos(dp)*radcon

	end function calculateAngle

	! Numerical Recipes material starts here
	! (C) Copr. 1986-92 Numerical Recipes Software

	! Gauss Jordan elimination to find matrix inverse
	subroutine gaussj(a,n,np)
	INTEGER n,np,NMAX
	DOUBLE PRECISION a(np,np)
	PARAMETER (NMAX=10)
	INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
	DOUBLE PRECISION big,dum,pivinv
	do 11 j=1,n
	  ipiv(j)=0
11	continue
	do 22 i=1,n
	  big=0.d0
	  do 13 j=1,n
	    if(ipiv(j).ne.1)then
	      do 12 k=1,n
	        if (ipiv(k).eq.0) then
	          if (abs(a(j,k)).ge.big)then
	            big=abs(a(j,k))
	            irow=j
	            icol=k
	          endif
	        else if (ipiv(k).gt.1) then
	          stop "Singular matrix encountered"
	        endif
12	      continue
	    endif
13	  continue
	  ipiv(icol)=ipiv(icol)+1
	  if (irow.ne.icol) then
	    do 14 l=1,n
	      dum=a(irow,l)
	      a(irow,l)=a(icol,l)
	      a(icol,l)=dum
14	    continue
	  endif
	  indxr(i)=irow
	  indxc(i)=icol
	  if (a(icol,icol).eq.0.d0) stop "Singular matrix encountered"
	  pivinv=1.d0/a(icol,icol)
	  a(icol,icol)=1.d0
	  do 16 l=1,n
	    a(icol,l)=a(icol,l)*pivinv
16	  continue
	  do 21 ll=1,n
	    if(ll.ne.icol)then
	      dum=a(ll,icol)
	      a(ll,icol)=0.d0
	      do 18 l=1,n
	        a(ll,l)=a(ll,l)-a(icol,l)*dum
18	      continue
	    endif
21	  continue
22	continue
	do 24 l=n,1,-1
	  if(indxr(l).ne.indxc(l))then
	    do 23 k=1,n
	      dum=a(k,indxr(l))
	      a(k,indxr(l))=a(k,indxc(l))
	      a(k,indxc(l))=dum
23	    continue
	  endif
24	continue
	return
	end subroutine

	end module utility
