!	** gluchb **
!	Outputs the 'best' species2 involved in hydrogen bonding with species1/mol1/atom[s]
!	Also 

	program gluchb
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,resfile
	character*20 :: temparg
	integer :: n,nframes,success,nargs,m,i,unique,zero
	integer :: framestodo, nsites, t(9)
	integer :: species2(3), sp2 = 2
	integer :: iargc
	integer, allocatable :: sites(:,:), id(:), idtemp(:)
	real*8 :: a(3),b(3),c(3),temp(3),v1(3),v2(3),v1mag,v2mag,roo,angle,dp,maxoo,minang,totaverage,minoo
	real*8, allocatable :: best(:)
	real*8 :: calcangle, delta

	minoo = 2.0

	nargs = iargc()
	if (nargs.LT.9) stop "Usage : pickhb <HISTORYfile> <OUTPUTfile> <nframes> <maxoo> <minang> <sp2> <O> <H1> <H2>"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temparg); read(temparg,"(I5)") framestodo
	call getarg(4,temparg); read(temparg,"(F12.4)") maxoo
	call getarg(5,temparg); read(temparg,"(F12.4)") minang
	call getarg(6,temparg); read(temparg,"(I4)") sp2
	call getarg(7,temparg); read(temparg,"(I4)") species2(1)
	call getarg(8,temparg); read(temparg,"(I4)") species2(2)
	call getarg(9,temparg); read(temparg,"(I4)") species2(3)
	
	write(0,*) "O...site distance min/max :",minoo,maxoo
	write(0,*) "Minimum O-H...site angle  :",minang

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	nsites = 6

	allocate(sites(nsites,2))
	allocate(id(nsites))
	allocate(idtemp(nsites))
	allocate(best(nsites))

	sites(:,:) = reshape( (/ 6,8,9,10,11,12,0,20,21,22,23,24 /) , (/ nsites,2 /) )

	write(0,"(A,4I3)") "'Water' species, O, H, H :",sp2,species2
	write(0,*) "OH site specification (O, H):"
	do m=1,nsites
	  write(0,*) (sites(m,n),n=1,2)
	end do

	open(unit=9,file="pickhb.dat",form="formatted",status="replace")
	write(9,"(A)") "   Frame   dO8   dO9  dO10  dO11  dO12  dO13  "

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	best = 100.0
	id = 0

	i = s_start(sp2)-1
	do n=1,s_nmols(sp2)
	  c(1)=xpos(i+species2(1))		! Water oxygen
	  c(2)=ypos(i+species2(1))
	  c(3)=zpos(i+species2(1))
	  do m=1,nsites
	    a(1)=xpos(sites(m,1))	! Glucose Oxygen
	    a(2)=ypos(sites(m,1))
	    a(3)=zpos(sites(m,1))
	    ! First, search for H-Bonds where glucose OH is the donor (glucO-glucH...waterO)
	    if (sites(m,2).NE.0) then
	      b(1)=xpos(sites(m,2))	! Glucose Hydrogen
	      b(2)=ypos(sites(m,2))
	      b(3)=zpos(sites(m,2))
	      ! gO-wO distance
	      call pbc(a(1),a(2),a(3),c(1),c(2),c(3),temp(1),temp(2),temp(3))
	      roo=sqrt( (temp(1)-c(1))**2 + (temp(2)-c(2))**2 + (temp(3)-c(3))**2 )
	      ! gO-gH-wO angle
              angle = calcangle(a,b,c)
	      if ((roo.LT.maxoo).AND.(angle.GT.minang)) then
		! delta = (roo * roo)
		delta = (roo*roo) + ((180.0-angle)/100.0)**2
		if (delta.LT.best(m)) then
		  best(m) = delta
		  id(m) = n
		end if
	      end if
	    end if
	    ! Now for H-bonds where the glucose O is the acceptor (glucO...waterH-waterO)
	    if (species2(2).NE.0) then
	      b(1)=xpos(i+species2(2))	! Water Hydrogen 1
	      b(2)=ypos(i+species2(2))
	      b(3)=zpos(i+species2(2))
	      ! gO-wO distance (calc again to be safe)
	      call pbc(a(1),a(2),a(3),c(1),c(2),c(3),temp(1),temp(2),temp(3))
	      roo=sqrt( (temp(1)-c(1))**2 + (temp(2)-c(2))**2 + (temp(3)-c(3))**2 )
	      ! gO-wH-wO angle
              angle = calcangle(a,b,c)
	      if ((roo.LT.maxoo).AND.(angle.GT.minang)) then
		! delta = (roo * roo)
		delta = (roo*roo) + ((180.0-angle)/100.0)**2
		if (delta.LT.best(m)) then
		  best(m) = delta
		  id(m) = n
		end if
	      end if
	      b(1)=xpos(i+species2(3))	! Water Hydrogen 2
	      b(2)=ypos(i+species2(3))
	      b(3)=zpos(i+species2(3))
	      ! O-O Dist is still the same, so just need angle.
              angle = calcangle(a,b,c)
	      if ((roo.LT.maxoo).AND.(angle.GT.minang)) then
		! delta = (roo * roo)
		delta = (roo*roo) + ((180.0-angle)/100.0)**2
		if (delta.LT.best(m)) then
		  best(m) = delta
		  id(m) = n
		end if
	      end if
	    end if
	  end do
	  i = i+s_natoms(sp2)
	end do

	zero = 0
	do n=1,nsites
	  if (id(n).EQ.0) zero = zero + 1
	end do

	idtemp = id
	unique = nsites - zero 
	do n=1,nsites-1
	  if (idtemp(n).NE.0) then
	    do m=n+1,nsites
	      if (idtemp(n).EQ.idtemp(m))  then
		unique = unique - 1
		idtemp(m) = 0
	      end if
	    end do
	  end if
	end do

	write(9,"(I5,2x,10I6)") nframes,(id(m),m=1,nsites),unique,zero

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

	close(9)

	write(0,*) "Finished."
999	close(14)
	end program gluchb

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
        calcangle = safeAngle(dp)
	end function calcangle
