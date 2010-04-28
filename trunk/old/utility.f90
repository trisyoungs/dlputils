	! Useful routines for many other programs!

	module utility
	  ! Centre-of-mass arrays
	  real*8, allocatable :: comx(:,:),comy(:,:),comz(:,:)
	contains

	subroutine pbc(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	use dlprw; implicit none
	! Performs minimum image convention.
	! xyz1: point to consider; xyz2: reference point; xyz3: minimum image coords of xyz1 (result)
	real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3
	! Performs minimum image convention.....
	x3=x1 - cell(1)*NINT((x1-x2)/cell(1))
	y3=y1 - cell(5)*NINT((y1-y2)/cell(5))
	z3=z1 - cell(9)*NINT((z1-z2)/cell(9))
	end subroutine pbc

	subroutine calc_com
	use dlprw; implicit none
	! Calculate the centre of mass for each molecule of all species
	real*8 :: massnorm,tx,ty,tz,x0,y0,z0
	integer :: s,sp,n,m,o,p
	do n=1,s_nmols(sp)
	  comx(sp,n)=0.0
	  comy(sp,n)=0.0
	  comz(sp,n)=0.0
	end do
	do sp=1,nspecies
	  if (s_natoms(sp).GT.1) THEN    ! Do centre-of-mass calculation
	    o=s_start(sp)
	    do n=1,s_nmols(sp)
	      ! We will use [0,0,0] as the reference point....
	      x0=0.0
	      y0=0.0
	      z0=0.0
	      massnorm=0.0
	      do m=1,s_natoms(sp)
	        CALL PBC(x(o+m-1),y(o+m-1),z(o+m-1),comx(sp,n),comy(sp,n),comz(sp,n),tx,ty,tz,cell)
	        x0=x0 + (tx-comx(sp,n))*mass(o+m-1)
	        y0=y0 + (ty-comy(sp,n))*mass(o+m-1)
	        z0=z0 + (tz-comz(sp,n))*mass(o+m-1)
	        massnorm=massnorm+mass(o+m-1)
	      end do
	      ! Normalise...
	      comx(sp,n)=comx(sp,n) + x0/massnorm
	      comy(sp,n)=comy(sp,n) + y0/massnorm
	      comz(sp,n)=comz(sp,n) + z0/massnorm 
	      o=o+s_natoms(sp)
	    end do
	  end if
	end do
	end subroutine calc_com

	end module utility
