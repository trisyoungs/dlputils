	! Useful routines for many other programs!

	module utility
	  ! Centre-of-mass arrays
	  real*8, allocatable :: comx(:,:),comy(:,:),comz(:,:)
	  ! Axis Arrays
	  integer, allocatable :: aa(:,:)			! Axis atoms
	  logical, allocatable :: axisdefined(:), axisusecom(:)
	  real*8, allocatable :: axisx(:,:,:), axisy(:,:,:), axisz(:,:,:)   ! Axis point coords
	  real*8, allocatable :: axisox(:,:), axisoy(:,:), axisoz(:,:)	! Axis origin
	  logical :: com_alloc = .false., axis_alloc = .false.
	  ! Reciprocal cell
	  real*8 :: rcell(9), celldet
	  ! Inverse cell
	  real*8 :: icell(9)
	contains

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

	subroutine calc_rcell()
	use dlprw; implicit none
	real*8 :: vol
	real*8, parameter :: pi = 3.14159265358979d0
	integer :: n
	if (imcon.gt.3) stop "calc_rcell can't handle this image convention."
!	rcell(1)=cell(5)*cell(9)-cell(6)*cell(8)
!	rcell(2)=cell(6)*cell(7)-cell(4)*cell(9)
!	rcell(3)=cell(4)*cell(8)-cell(5)*cell(7)
!	rcell(4)=cell(3)*cell(8)-cell(2)*cell(9)
!	rcell(5)=cell(1)*cell(9)-cell(3)*cell(7)
!	rcell(6)=cell(2)*cell(7)-cell(1)*cell(8)
!	rcell(7)=cell(2)*cell(6)-cell(3)*cell(5)
!	rcell(8)=cell(3)*cell(4)-cell(1)*cell(6)
!	rcell(9)=cell(1)*cell(5)-cell(2)*cell(4)
!	vol = abs( cell(1)*rcell(1) + cell(2)*rcell(2) + cell(3)*rcell(3))
!	do n=1,9
!	  rcell(n)=2.0d0*pi*rcell(n)/vol
!	end do

	! Calculate the reciprocal cell (cubic, orthorhombic, or parallelepiped)
	! Take cross products of original cell vectors to get reciprocal cell vectors
	! R(1,4,7) = C(Y) * C(Z)
	! R(2,5,8) = C(X) * C(Z)
	! R(3,6,9) = C(X) * C(Y)
	rcell(1) = cell(5) * cell(9) - cell(8) * cell(6)
	rcell(4) = cell(8) * cell(3) - cell(2) * cell(9)
	rcell(7) = cell(2) * cell(6) - cell(5) * cell(3)
	rcell(2) = cell(4) * cell(9) - cell(7) * cell(6)
	rcell(5) = cell(7) * cell(3) - cell(1) * cell(9)
	rcell(8) = cell(1) * cell(6) - cell(4) * cell(3)
	rcell(3) = cell(4) * cell(8) - cell(7) * cell(5)
	rcell(6) = cell(7) * cell(2) - cell(1) * cell(8)
	rcell(9) = cell(1) * cell(5) - cell(4) * cell(2)
	vol = dabs(cell(1)*rcell(1) + cell(5)*rcell(5) + cell(9)*rcell(9))
	rcell = 2.0d0*pi*rcell / vol
	write(0,"(3F12.4)") rcell
	end subroutine calc_rcell

	subroutine calc_icell()
	use dlprw; implicit none
	icell = cell
	call gaussj(icell,3,3)
	has_icell = .true.
	end subroutine calc_icell

	subroutine genaxis
	use dlprw;
	implicit none
	! Subroutine to calculate individual molecular axes of species molecules
	real*8 :: vx,vy,vz,tx,ty,tz
	real*8 :: xmag, ymag, dp
	integer :: sp,o,n
	if (.NOT.axis_alloc) call alloc_axis
	do sp=1,nspecies
	  ! Calculate local axes if we can
	  if (axisdefined(sp).and.(s_natoms(sp).gt.2).and.(aa(sp,1).ne.aa(sp,2))) then
	    o=s_start(sp)
	    do n=1,s_nmols(sp)
	      ! 1) Determine the X axis of the molecule
	      !    -- Get the vector between atoms 1 and 2...
	      call pbc(xpos(o+aa(sp,2)-1),ypos(o+aa(sp,2)-1),zpos(o+aa(sp,2)-1), &
	      & xpos(o+aa(sp,1)-1),ypos(o+aa(sp,1)-1),zpos(o+aa(sp,1)-1),tx,ty,tz)
	      call getvector(tx,ty,tz,xpos(o+aa(sp,1)-1),ypos(o+aa(sp,1)-1),zpos(o+aa(sp,1)-1),vx,vy,vz)
	      xmag=SQRT(vx**2 + vy**2 + vz**2)     ! Magnitude of vector
	      axisx(sp,n,1)=vx		! Store the un-normalised components of the x-axis vector
	      axisx(sp,n,2)=vy
	      axisx(sp,n,3)=vz
	      ! 2) Set the origin of the axis system - halfway along the x-axis vector
	      axisox(sp,n)=xpos(o+aa(sp,1)-1) + 0.5*vx		! N.B. *NOT* the normalised axis, but the actual distance i->j
	      axisoy(sp,n)=ypos(o+aa(sp,1)-1) + 0.5*vy
	      axisoz(sp,n)=zpos(o+aa(sp,1)-1) + 0.5*vz
	      ! 3) Determine the Y axis of the molecule
	      !    -- Get the vector between midway along the x-axis vector and midway between
	      !	the positions of aa 3 and 4.
	      call pbc(xpos(o+aa(sp,3)-1),ypos(o+aa(sp,3)-1),zpos(o+aa(sp,3)-1), &
	      & axisox(sp,n),axisoy(sp,n),axisoz(sp,n),tx,ty,tz)
	      call pbc(xpos(o+aa(sp,4)-1),ypos(o+aa(sp,4)-1),zpos(o+aa(sp,4)-1), &
	      & axisox(sp,n),axisoy(sp,n),axisoz(sp,n),vx,vy,vz)	! Use v_xyz temporarily
	      tx = 0.5*(tx+vx)
	      ty = 0.5*(ty+vy)
	      tz = 0.5*(tz+vz)
	      call getvector(axisox(sp,n),axisoy(sp,n),axisoz(sp,n),tx,ty,tz,vx,vy,vz)
	      axisy(sp,n,1)=vx   ! Store the un-normalised components of the y-axis vector
	      axisy(sp,n,2)=vy
	      axisy(sp,n,3)=vz
	      ! 3a) Orthogonalise the y vector w.r.t. the x vector via Gram-Schmidt
	      dp = axisx(sp,n,1)*axisy(sp,n,1) + axisx(sp,n,2)*axisy(sp,n,2) + axisx(sp,n,3)*axisy(sp,n,3)
	      tx = axisy(sp,n,1) - (dp / xmag**2)*axisx(sp,n,1)
	      ty = axisy(sp,n,2) - (dp / xmag**2)*axisx(sp,n,2)
	      tz = axisy(sp,n,3) - (dp / xmag**2)*axisx(sp,n,3)
	      axisy(sp,n,1) = tx
	      axisy(sp,n,2) = ty
	      axisy(sp,n,3) = tz
	      ymag=SQRT(tx**2 + ty**2 + tz**2)     ! Magnitude of vector
	      ! 3b) Normalise the y and x vectors
	      axisx(sp,n,:) = axisx(sp,n,:) / xmag
	      axisy(sp,n,:) = axisy(sp,n,:) / ymag
	      ! 4) Calculate the z axis from the x and y axes....
	      axisz(sp,n,1)=axisx(sp,n,2)*axisy(sp,n,3) - axisx(sp,n,3)*axisy(sp,n,2)
	      axisz(sp,n,2)=axisx(sp,n,3)*axisy(sp,n,1) - axisx(sp,n,1)*axisy(sp,n,3)
	      axisz(sp,n,3)=axisx(sp,n,1)*axisy(sp,n,2) - axisx(sp,n,2)*axisy(sp,n,1)
	      ! Increase the atom count for the next molecule....
	      o=o+s_natoms(sp)
	      !write(88,*) "O: ",axisox(sp,n),axisoy(sp,n),axisoz(sp,n)
	      !write(88,*) "X: ",axisx(sp,n,1),axisx(sp,n,2),axisx(sp,n,3)
	      !write(88,*) "Y: ",axisy(sp,n,1),axisy(sp,n,2),axisy(sp,n,3)
	      !write(88,*) "Z: ",axisz(sp,n,1),axisz(sp,n,2),axisz(sp,n,3)
	    end do
	  else
	    ! No axis definition, or not enough atoms, or two identical x-axis IDs supplied
	    ! Set the axes to be in standard orientation
	    if ((.not.axisdefined(sp)).and.(axisusecom(sp))) call calc_spcom(sp)
	    o=s_start(sp)
	    do n=1,s_nmols(sp)
	      axisx(sp,n,:) = (/ 1.0,0.0,0.0 /)
	      axisy(sp,n,:) = (/ 0.0,1.0,0.0 /)
	      axisz(sp,n,:) = (/ 0.0,0.0,1.0 /)
	      if (axisusecom(sp)) then
	        axisox(sp,n)=comx(sp,n)
	        axisoy(sp,n)=comy(sp,n)
	        axisoz(sp,n)=comz(sp,n)
	      else
	        axisox(sp,n)=xpos(o+aa(sp,1)-1)
	        axisoy(sp,n)=ypos(o+aa(sp,1)-1)
	        axisoz(sp,n)=zpos(o+aa(sp,1)-1)
	      end if
	      o=o+s_natoms(sp)
	    end do
	  end if
	end do
	end subroutine genaxis

	subroutine alloc_axis()
	use dlprw; implicit none
	integer :: status
	axis_alloc = .true.
	allocate(axisx(nspecies,maxmols,3),stat=status); if (status.GT.0) stop "Allocation error for axisx()"
	allocate(axisy(nspecies,maxmols,3),stat=status); if (status.GT.0) stop "Allocation error for axisy()"
	allocate(axisz(nspecies,maxmols,3),stat=status); if (status.GT.0) stop "Allocation error for axisz()"
	allocate(axisox(nspecies,maxmols),stat=status); if (status.GT.0) stop "Allocation error for axisox()"
	allocate(axisoy(nspecies,maxmols),stat=status); if (status.GT.0) stop "Allocation error for axisoy()"
	allocate(axisoz(nspecies,maxmols),stat=status); if (status.GT.0) stop "Allocation error for axisoz()"
	allocate(aa(nspecies,4))
	allocate(axisdefined(nspecies))
	allocate(axisusecom(nspecies))
	axisdefined = .false.
	axisusecom = .true.
	aa = 0
	end subroutine alloc_axis

	subroutine getvector(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	implicit none
	real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3
	! xyz1,xyz2: coordinates of two atoms; xyz3: vector between them (result)
	x3=x1-x2
	y3=y1-y2
	z3=z1-z2
	end subroutine getvector

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

	real*8 function volume(cell)
	implicit none
	real*8, intent(in) :: cell(9)
	! Hard-coded calculation of determinant of 3x3 matrix
	volume = cell(1) * (cell(5)*cell(9) - cell(8)*cell(6));
	volume = volume - cell(2) * (cell(4)*cell(9) - cell(7)*cell(6));
	volume = volume + cell(3) * (cell(4)*cell(8) - cell(7)*cell(5));
	end function volume

	! Numerical Recipes material starts here

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
