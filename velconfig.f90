	! Re-velocitise a DL_POLY config with a Gaussian distribution at a certain T

	program velconfig
	use dlprw; use utility; use dlpfield
	implicit none
	character*80 :: configfile, temp, fieldfile
	real*8 :: t, rr, rr1, rr2
	real*8, parameter :: pi = 3.14159265358979d0
	integer :: nargs,n,i,m,success,seed,sp
	real*8 :: duni, ran2
	real*8 :: sigma,cmx,cmy,cmz,cmvx,cmvy,cmvz,sysmas, degfre
	real*8 :: amx,amy,amz,det,scale,rsq,wxx,wyy,wzz,sumke,roti(9),rotinv(9)

	nargs = iargc()
	if (nargs.ne.4) stop "Usage: velconfig <CONFIG file> <FIELD file> <temperature, in K> <random seed, or 0 for DL_POLY values>"
	call getarg(1,configfile)
	call getarg(2,fieldfile)
	call getarg(3,temp)
	read(temp,"(f10.5)") t
	call getarg(4,temp)
	read(temp,"(i10)") seed

	! Read FIELD file
	if (.not.readfield(fieldfile)) stop "Couldn't read FIELD file."
	natms = f_totalatoms
	call alloc_xyz()
	
	! Open the config file
	call readconfig(configfile)
	write(0,*) "Number of atoms in CONFIG file: ", natms
	keytrj = 1

	! Initialise random number generator and assign random velocities
	if (seed.eq.0) then
	  rr = duni()

	  do i=1,2*(natms/2),2
	  
	   rr1=sqrt(-2.d0*log(duni()))
	   rr2=2.d0*pi*duni()
	   cfgvel(i,1)=rr1*cos(rr2)
	   cfgvel(i+1,1)=rr1*sin(rr2)

	   rr1=sqrt(-2.d0*log(duni()))
	   rr2=2.d0*pi*duni()
	   cfgvel(i,2)=rr1*cos(rr2)
	   cfgvel(i+1,2)=rr1*sin(rr2)

	   rr1=sqrt(-2.d0*log(duni()))
	   rr2=2.d0*pi*duni()
	   cfgvel(i,3)=rr1*cos(rr2)
	   cfgvel(i+1,3)=rr1*sin(rr2)
	  
	  enddo
	  if(mod(natms,2).ne.0)then
	  
	    rr1=sqrt(-2.d0*log(duni()))
	    rr2=2.d0*pi*duni()
	    cfgvel(natms,1)=rr1*cos(rr2)
	    cfgvel(natms,2)=rr1*sin(rr2)
	    rr1=sqrt(-2.d0*log(duni()))
	    rr2=2.d0*pi*duni()
	    cfgvel(natms,3)=rr1*cos(rr2)
	  
	  endif
	else
	  rr = ran2(seed)

	  do i=1,2*(natms/2),2
	  
	   rr1=sqrt(-2.d0*log(ran2(seed)))
	   rr2=2.d0*pi*ran2(seed)
	   cfgvel(i,1)=rr1*cos(rr2)
	   cfgvel(i+1,1)=rr1*sin(rr2)

	   rr1=sqrt(-2.d0*log(ran2(seed)))
	   rr2=2.d0*pi*ran2(seed)
	   cfgvel(i,2)=rr1*cos(rr2)
	   cfgvel(i+1,2)=rr1*sin(rr2)

	   rr1=sqrt(-2.d0*log(ran2(seed)))
	   rr2=2.d0*pi*ran2(seed)
	   cfgvel(i,3)=rr1*cos(rr2)
	   cfgvel(i+1,3)=rr1*sin(rr2)
	  
	  enddo
	  if(mod(natms,2).ne.0)then
	  
	    rr1=sqrt(-2.d0*log(ran2(seed)))
	    rr2=2.d0*pi*ran2(seed)
	    cfgvel(natms,1)=rr1*cos(rr2)
	    cfgvel(natms,2)=rr1*sin(rr2)
	    rr1=sqrt(-2.d0*log(ran2(seed)))
	    rr2=2.d0*pi*ran2(seed)
	    cfgvel(natms,3)=rr1*cos(rr2)
	  end if  
	end if

	! dl_poly subroutine for scaling the velocity arrays to the
	! desired temperature
	! zeroes angular momentum in non-periodic system.
      
	degfre = f_totalatoms * 3
	do n=1,f_nspecies
	  degfre = degfre - f_nmols(n)*f_nconstraints(n)
	end do
	degfre = degfre - 3
	if (imcon.eq.0) degfre = degfre - 3
	write(0,*) "Degrees of freedom: ",degfre
	
	sigma = t*8.31451115d-1*degfre*0.5d0

	! calculate centre of mass position and motion of the system
	cmx=0.d0
	cmy=0.d0
	cmz=0.d0
	cmvx=0.d0
	cmvy=0.d0
	cmvz=0.d0
	sysmas=0.d0
	
	i=1
	do sp=1,f_nspecies
	do m=1,f_nmols(sp)
	do n=1,f_natoms(sp)
	  if(f_masses(sp,n).gt.1.d-6) then
	    cmx=cmx+cfgxyz(i,1)*f_masses(sp,n)
	    cmy=cmy+cfgxyz(i,2)*f_masses(sp,n)
	    cmz=cmz+cfgxyz(i,3)*f_masses(sp,n)
	    sysmas=sysmas+f_masses(sp,n)
	    cmvx=cmvx+cfgvel(i,1)*f_masses(sp,n)
	    cmvy=cmvy+cfgvel(i,2)*f_masses(sp,n)
	    cmvz=cmvz+cfgvel(i,3)*f_masses(sp,n)
	  endif
	  i = i + 1
	enddo
	enddo
	enddo
	
	cmx=cmx/sysmas
	cmy=cmy/sysmas
	cmz=cmz/sysmas
	
	cmvx=cmvx/sysmas
	cmvy=cmvy/sysmas
	cmvz=cmvz/sysmas
	
	! remove centre of mass motion  
	i=1
	do sp=1,f_nspecies
	do m=1,f_nmols(sp)
	do n=1,f_natoms(sp)
	  if(f_masses(sp,n).gt.1.d-6) then
	    cfgvel(i,1)=cfgvel(i,1)-cmvx
	    cfgvel(i,2)=cfgvel(i,2)-cmvy
	    cfgvel(i,3)=cfgvel(i,3)-cmvz
	  else
	    cfgvel(i,1)=0.d0
	    cfgvel(i,2)=0.d0
	    cfgvel(i,3)=0.d0
	  endif
	  i = i + 1
	enddo
	enddo
	enddo
	
	! zero angular momentum about centre of mass - non-periodic system
	if(imcon.eq.0) then

	! move to centre of mass origin
	  do i=1,natms
	    cfgxyz(i,1)=cfgxyz(i,1)-cmx
	    cfgxyz(i,2)=cfgxyz(i,2)-cmy
	    cfgxyz(i,3)=cfgxyz(i,3)-cmz
	  enddo
	  
	! angular momentum accumulators
	  amx=0.d0
	  amy=0.d0
	  amz=0.d0

	! rotational inertia accumulators
	  do i=1,9
	    roti(i)=0.d0
	  enddo
	  
	  i=1
	  do sp=1,f_nspecies
	  do m=1,f_nmols(sp)
	  do n=1,f_natoms(sp)
	    amx=amx+f_masses(sp,n)*(cfgxyz(i,2)*cfgvel(i,3)-cfgxyz(i,3)*cfgvel(i,2))
	    amy=amy+f_masses(sp,n)*(cfgxyz(i,3)*cfgvel(i,1)-cfgxyz(i,1)*cfgvel(i,3))
	    amz=amz+f_masses(sp,n)*(cfgxyz(i,1)*cfgvel(i,2)-cfgxyz(i,2)*cfgvel(i,1))
	    
	    rsq=cfgxyz(i,1)**2+cfgxyz(i,2)**2+cfgxyz(i,3)**2
	    roti(1)=roti(1)+f_masses(sp,n)*(cfgxyz(i,1)*cfgxyz(i,1)-rsq)
	    roti(2)=roti(2)+f_masses(sp,n)* cfgxyz(i,1)*cfgxyz(i,2)
	    roti(3)=roti(3)+f_masses(sp,n)* cfgxyz(i,1)*cfgxyz(i,3)
	    roti(5)=roti(5)+f_masses(sp,n)*(cfgxyz(i,2)*cfgxyz(i,2)-rsq)
	    roti(6)=roti(6)+f_masses(sp,n)* cfgxyz(i,2)*cfgxyz(i,3)
	    roti(9)=roti(9)+f_masses(sp,n)*(cfgxyz(i,3)*cfgxyz(i,3)-rsq)
	    i = i + 1
	  enddo
	  enddo
	  enddo

	! complete rotational inertia matrix
	  roti(4)=roti(2)
	  roti(7)=roti(3)
	  roti(8)=roti(6)

	! invert rotational inertia matrix
	  call invert (roti,rotinv,det)

	! correction to angular velocity
	  wxx=rotinv(1)*amx+rotinv(2)*amy+rotinv(3)*amz
	  wyy=rotinv(4)*amx+rotinv(5)*amy+rotinv(6)*amz
	  wzz=rotinv(7)*amx+rotinv(8)*amy+rotinv(9)*amz

	! correction to linear velocity
	  i=1
	  do sp=1,f_nspecies
	  do m=1,f_nmols(sp)
	  do n=1,f_natoms(sp)
	    if(f_masses(sp,n).gt.1.d-6) then
	      cfgvel(i,1)=cfgvel(i,1)+(wyy*cfgxyz(i,3)-wzz*cfgxyz(i,2))
	      cfgvel(i,2)=cfgvel(i,2)+(wzz*cfgxyz(i,1)-wxx*cfgxyz(i,3))
	      cfgvel(i,3)=cfgvel(i,3)+(wxx*cfgxyz(i,2)-wyy*cfgxyz(i,1))
	    endif
	    i = i + 1
	  enddo
	  enddo
	  enddo

	! reset positions to original reference frame
	  do i=1,natms
	    cfgxyz(i,1)=cfgxyz(i,1)+cmx
	    cfgxyz(i,2)=cfgxyz(i,2)+cmy
	    cfgxyz(i,3)=cfgxyz(i,3)+cmz
	  enddo
	  
	endif

	! calculate temperature 
	sumke=0.d0
	i=1
	do sp=1,f_nspecies
	do m=1,f_nmols(sp)
	do n=1,f_natoms(sp)
	  sumke=sumke+f_masses(sp,n)*(cfgvel(i,1)**2+cfgvel(i,2)**2+cfgvel(i,3)**2)
	  i = i + 1
	enddo
	enddo
	enddo

	sumke=0.5d0*sumke

	! apply temperature scaling
	scale=1.d0
	if(sumke.gt.1.d-6)scale=sqrt(sigma/sumke)
	do i=1,natms
	  cfgvel(i,1)=scale*cfgvel(i,1)
	  cfgvel(i,2)=scale*cfgvel(i,2)
	  cfgvel(i,3)=scale*cfgvel(i,3)
	enddo

	call writeconfig(configfile)

	stop "Finished."

	end program velconfig


	function duni()

	!*********************************************************************
	!     
	!     dl_poly random number generator based on the universal
	!     random number generator of marsaglia, zaman and tsang
	!     (stats and prob. lett. 8 (1990) 35-39.) it must be
	!     called once to initialise parameters u,c,cd,cm
	!     
	!     copyright daresbury laboratory 1992
	!     author -  w.smith	   july 1992
	!     
	!     wl
	!     2008/12/23 10:29:12
	!     1.9
	!     Exp
	!     
	!*********************************************************************

	implicit none

	logical new
	integer ir,jr,i,j,k,l,m,ii,jj
	real(4) s,t,u,c,cd,cm,uni
	real(8) duni
	dimension u(97)
	save u,c,cd,cm,uni,ir,jr,new
	data new/.true./

	if(new)then

	!     initial values of i,j,k must be in range 1 to 178 (not all 1)
	!     initial value of l must be in range 0 to 168.

	  i=12
	  j=34
	  k=56
	  l=78
     
	  ir=97
	  jr=33
	  new=.false.

	  do 200 ii=1,97
	    s=0.0
	    t=0.5
	    do 100 jj=1,24
	      m=mod(mod(i*j,179)*k,179)
	      i=j
	      j=k
	      k=m
	      l=mod(53*l+1,169)
	      if(mod(l*m,64).ge.32)s=s+t
	      t=0.5*t
100	continue
	    u(ii)=s
200	continue
	  c =  362436.0/16777216.0
	  cd= 7654321.0/16777216.0
	  cm=16777213.0/16777216.0
	else

	! calculate random number
	  uni=u(ir)-u(jr)
	  if(uni.lt.0.0)uni=uni+1.0
	  u(ir)=uni
	  ir=ir-1
	  if(ir.eq.0)ir=97
	  jr=jr-1
	  if(jr.eq.0)jr=97
	  c=c-cd
	  if(c.lt.0.0)c=c+cm
	  uni=uni-c
	  if(uni.lt.0.0)uni=uni+1.0
	  duni=dble(uni)
	endif
	return
	end function duni

	
	real*8 function ran2(idum)
	implicit none
	integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
	real*8 :: AM,EPS,RNMX
	parameter (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1, &
	& IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
	& NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
	! 'ran2' random number generator.
	!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.
	integer :: idum2,j,k,iv(NTAB),iy
	save iv,iy,idum2
	data idum2/123456789/, iv/NTAB*0/, iy/0/
	if (idum.le.0) then
	  idum=max(-idum,1)
	  idum2=idum
	  do 11 j=NTAB+8,1,-1
	    k=idum/IQ1
	    idum=IA1*(idum-k*IQ1)-k*IR1
	    if (idum.lt.0) idum=idum+IM1
	    if (j.le.NTAB) iv(j)=idum
11	continue
	  iy=iv(1)
	endif
	k=idum/IQ1
	idum=IA1*(idum-k*IQ1)-k*IR1
	if (idum.lt.0) idum=idum+IM1
	k=idum2/IQ2
	idum2=IA2*(idum2-k*IQ2)-k*IR2
	if (idum2.lt.0) idum2=idum2+IM2
	j=1+iy/NDIV
	iy=iv(j)-idum2
	iv(j)=idum
	if(iy.lt.1)iy=iy+IMM1
	ran2=min(AM*iy,RNMX)
	return
	end function ran2
	
	subroutine invert(a,b,d)
	! dl_poly subroutine to invert a 3 * 3 matrix using cofactors
	! copyright - daresbury laboratory 1992
	! author    - w. smith	 april 1992
	implicit none

	real(8) a,b,d,r

	dimension a(9),b(9)

	! calculate adjoint matrix
	b(1)=a(5)*a(9)-a(6)*a(8)
	b(2)=a(3)*a(8)-a(2)*a(9)
	b(3)=a(2)*a(6)-a(3)*a(5)
	b(4)=a(6)*a(7)-a(4)*a(9)
	b(5)=a(1)*a(9)-a(3)*a(7)
	b(6)=a(3)*a(4)-a(1)*a(6)
	b(7)=a(4)*a(8)-a(5)*a(7)
	b(8)=a(2)*a(7)-a(1)*a(8)
	b(9)=a(1)*a(5)-a(2)*a(4)

	! calculate determinant
	d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
	r=0.d0
	if(abs(d).gt.0.d0)r=1.d0/d

	! complete inverse matrix
	b(1)=r*b(1)
	b(2)=r*b(2)
	b(3)=r*b(3)
	b(4)=r*b(4)
	b(5)=r*b(5)
	b(6)=r*b(6)
	b(7)=r*b(7)
	b(8)=r*b(8)
	b(9)=r*b(9)

	return
	end subroutine invert
