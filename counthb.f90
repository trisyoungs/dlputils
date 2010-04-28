!	** counthb **
!	Counts the number of h-bonds between specified species per frame.

	program counthb
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,resfile
	character*20 :: temparg
	integer :: n,nframes,success,nargs,m,i
	integer :: framestodo, nspec, t(9)
	integer :: glucoh(6,2), water(3), watersp = 2
	integer :: dgrid, agrid, nacceptor, ndonor, abin, dbin
	integer :: iargc
	real*8 :: a(3),b(3),c(3),temp(3),v1(3),v2(3),v1mag,v2mag,roo,angle,dp,maxoo,minang,totaverage,minoo
	integer, allocatable :: results(:,:,:), totals(:)
	real*8, allocatable :: frameaverages(:), averages(:,:)
	real*8, allocatable :: donordens(:,:), acceptordens(:,:)
	integer, parameter :: DONOR = 1, ACCEPTOR = 2
	logical :: printgeom = .FALSE.
	real*8 :: calcangle, ddelta = 0.01, adelta = 1.0

	minoo = 2.0

	nargs = iargc()
	if (nargs.LT.9) stop "Usage : counthb <DLP HISTORYfile> <DLP OUTPUTfile> <nframes> <maxoo> <minang> <wsp> <wO> <wH1> <wH2> [print]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temparg); read(temparg,"(I5)") framestodo
	call getarg(4,temparg); read(temparg,"(F12.4)") maxoo
	call getarg(5,temparg); read(temparg,"(F12.4)") minang
	call getarg(6,temparg); read(temparg,"(I4)") watersp
	call getarg(7,temparg); read(temparg,"(I4)") water(1)
	call getarg(8,temparg); read(temparg,"(I4)") water(2)
	call getarg(9,temparg); read(temparg,"(I4)") water(3)
	
	if (nargs.EQ.10) then
	  call getarg(10,temparg)
	  if (temparg.EQ."print") printgeom = .TRUE.
	end if

	write(0,*) "O-O distance min/max :",minoo,maxoo
	write(0,*) "Minimum O-H-O angle  :",minang

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	allocate(results(framestodo,6,2))
	allocate(totals(framestodo))
	allocate(frameaverages(framestodo))
	allocate(averages(framestodo,2))
	glucoh(:,:) = reshape( (/ 6,8,9,10,11,12,0,20,21,22,23,24 /) , (/ 6,2 /) )
	results(:,:,:) = reshape( (/ (0,n=1,2*6*framestodo) /) , (/ framestodo,6,2 /) )
	totals(:) = (/ (0.0,n=1,framestodo) /)
	totaverage = 0.0
	frameaverages = 0.0
	averages = 0.0

	write(0,"(A,4I3)") "Water species, O, H, H :",watersp,water
	write(0,*) "Glucose OH specification:"
	do m=1,6
	  write(0,*) (glucoh(m,n),n=1,2)
	end do

	! Prepare if printing is required
	if (printgeom) then
	  open(unit=14,file="counthb.geom",form="formatted",status="replace")
	  dgrid = int((maxoo-minoo) / ddelta)+1
	  agrid = int((180.0-minang) / adelta)+1
	  write(0,*) "Grid sizes (roo, angle) are ",dgrid,agrid
	  allocate(donordens(1:dgrid,1:agrid))
	  allocate(acceptordens(1:dgrid,1:agrid))
	  donordens = 0
	  acceptordens = 0
	  ndonor = 0
	  nacceptor = 0
	end if

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.EQ.-1) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	i = s_start(watersp)-1
	do n=1,s_nmols(watersp)
	  c(1)=xpos(i+water(1))		! Water oxygen
	  c(2)=ypos(i+water(1))
	  c(3)=zpos(i+water(1))
	  do m=1,6
	    a(1)=xpos(glucoh(m,1))	! Glucose Oxygen
	    a(2)=ypos(glucoh(m,1))
	    a(3)=zpos(glucoh(m,1))
	    ! First, search for H-Bonds where glucose OH is the donor (glucO-glucH...waterO)
	    if (glucoh(m,2).NE.0) then
	      b(1)=xpos(glucoh(m,2))	! Glucose Hydrogen
	      b(2)=ypos(glucoh(m,2))
	      b(3)=zpos(glucoh(m,2))
	      ! gO-wO distance
	      call pbc(a(1),a(2),a(3),c(1),c(2),c(3),temp(1),temp(2),temp(3))
	      roo=sqrt( (temp(1)-c(1))**2 + (temp(2)-c(2))**2 + (temp(3)-c(3))**2 )
	      ! gO-gH-wO angle
              angle = calcangle(a,b,c)
	      if ((roo.LT.maxoo).AND.(angle.GT.minang)) then
		results(nframes,m,DONOR) = results(nframes,m,DONOR) + 1.0
		if (printgeom) then
		  write(14,"(F7.5,2x,F8.3)") roo,angle
		  ndonor = ndonor + 1
		  dbin = nint((roo-minoo)/ddelta); abin = nint((angle-minang)/adelta)
		  donordens(dbin,abin) = donordens(dbin,abin) + 1
		end if
	      end if
	    end if
	    ! Now for H-bonds where the glucose O is the acceptor (glucO...waterH-waterO)
	    if (water(2).NE.0) then
	      b(1)=xpos(i+water(2))	! Water Hydrogen 1
	      b(2)=ypos(i+water(2))
	      b(3)=zpos(i+water(2))
	      ! gO-wO distance (calc again to be safe)
	      call pbc(a(1),a(2),a(3),c(1),c(2),c(3),temp(1),temp(2),temp(3))
	      roo=sqrt( (temp(1)-c(1))**2 + (temp(2)-c(2))**2 + (temp(3)-c(3))**2 )
	      ! gO-wH-wO angle
              angle = calcangle(a,b,c)
	      if ((roo.LT.maxoo).AND.(angle.GT.minang)) then
		results(nframes,m,ACCEPTOR) = results(nframes,m,ACCEPTOR) + 1.0
		if (printgeom) then
		  write(14,"(F7.5,12x,F8.3)") roo,angle
		  nacceptor = nacceptor + 1
		  dbin = nint((roo-minoo)/ddelta); abin = nint((angle-minang)/adelta)
		  acceptordens(dbin,abin) = acceptordens(dbin,abin) + 1
		end if
	      end if
	      b(1)=xpos(i+water(3))	! Water Hydrogen 2
	      b(2)=ypos(i+water(3))
	      b(3)=zpos(i+water(3))
	      ! O-O Dist is still the same, so just need angle.
              angle = calcangle(a,b,c)
	      if ((roo.LT.maxoo).AND.(angle.GT.minang)) then
		results(nframes,m,ACCEPTOR) = results(nframes,m,ACCEPTOR) + 1.0
		if (printgeom) then
		  write(14,"(F7.5,12x,F8.3)") roo,angle
		  nacceptor = nacceptor + 1
		  dbin = nint((roo-minoo)/ddelta); abin = nint((angle-minang)/adelta)
		  acceptordens(dbin,abin) = acceptordens(dbin,abin) + 1
		end if
	      end if
	    end if
	  end do
	  i = i+s_natoms(watersp)
	end do
	      
	! Create the totals for this frame
	do n=1,6
	  totals(nframes) = totals(nframes) + results(nframes,n,DONOR) + results(nframes,n,ACCEPTOR)
	  averages(n,DONOR) = averages(n,DONOR) + results(nframes,n,DONOR)
	  averages(n,ACCEPTOR) = averages(n,ACCEPTOR) + results(nframes,n,ACCEPTOR)
	end do
	totaverage = totaverage + real(totals(nframes))
	frameaverages(nframes) = totaverage / real(nframes)

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


	totaverage = totaverage / nframes
	do n=1,6
	  averages(n,DONOR) = averages(n,DONOR) / nframes
	  averages(n,ACCEPTOR) = averages(n,ACCEPTOR) / nframes
	end do
	
	open(unit=9,file="counthb.dat",form="formatted",status="replace")
	write(9,"(A)") "Frame   dO8   dO9  dO10  dO11  dO12  dO13  aO8   aO9  aO10  aO11  aO12  aO13  Total"
	do n=1,framestodo
	  write(9,"(I5,2x,13(I4,2X),F7.4)") n,(results(n,m,DONOR),m=1,6),(results(n,m,ACCEPTOR),m=1,6),totals(n),frameaverages(n)
	end do
	write(9,"(' Avg ',2x,13(F5.3,1X),F7.4)") (averages(n,DONOR),n=1,6),(averages(n,ACCEPTOR),n=1,6),totaverage
	close(9)

	if (printgeom) then
	  open(unit=9,file="counthb.ddens",form="formatted",status="replace")
	  donordens = donordens / ndonor
	  do n=1,dgrid
	    do m=1,agrid
	      write(9,"(F10.6)") donordens(n,m)
	    end do
	  end do
	  close(9)
	  open(unit=9,file="counthb.adens",form="formatted",status="replace")
	  acceptordens = acceptordens / nacceptor
	  do n=1,dgrid
	    do m=1,agrid
	      write(9,"(F10.6)") acceptordens(n,m)
	    end do
	  end do
	  close(9)
	end if

	write(0,*) "Finished."
999	close(14)
	end program counthb

	real*8 function calcangle(a,b,c)
	use utility; implicit none
	real*8,dimension(3) :: a,b,c,v1,v2,temp
	real*8 :: v1mag,v2mag,dp
        ! Angle between atoms a--b--c
	call pbc(a(1),a(2),a(3),b(1),b(2),b(3),temp(1),temp(2),temp(3))
	a = temp
	call pbc(c(1),c(2),c(3),b(1),b(2),b(3),temp(1),temp(2),temp(3))
	c = temp
	! Get the bond vectors
	v1 = b - a
	v2 = b - c
	v1mag = sqrt( v1(1)**2 + v1(2)**2 + v1(3)**2 )
	v2mag = sqrt( v2(1)**2 + v2(2)**2 + v2(3)**2 )
	v1 = v1 / v1mag
	v2 = v2 / v2mag
        ! Calculate dot product and angle...
        dp = (v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3))
        calcangle = acos(dp)* 57.29577951d0
	end function calcangle
