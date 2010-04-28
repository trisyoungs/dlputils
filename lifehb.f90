!	** lifehb **
!	Calculate the histograms of H-bond lifetimes over the simulation

	program lifehb
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,resfile
	character*20 :: temparg
	integer :: n,m,o,nframes,success,nargs,i,j,m1,m2,s
	integer :: framestodo, nsites
	integer :: sp1, species2(3), sp2, nhbonds
	integer :: iargc
	integer, allocatable :: sites(:,:), life(:,:)
	real*8 :: a(3),b(3),c(3),temp(3),v1(3),v2(3),v1mag,v2mag,roo,angle,dp,maxoo,minang,minoo
	real*8, allocatable :: hist(:)
	real*8 :: calcangle, framestep

	minoo = 2.0

	nargs = iargc()
	if (nargs.LT.11) stop "Usage : lifehb <DLP HISTORYfile> <DLP OUTPUTfile> <nframes> <maxoo> <minang> <sp2> <O> <H1> <H2> <framestep> <sp1> [<X> <H>]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temparg); read(temparg,"(I5)") framestodo
	call getarg(4,temparg); read(temparg,"(F12.4)") maxoo
	call getarg(5,temparg); read(temparg,"(F12.4)") minang
	call getarg(6,temparg); read(temparg,"(I4)") sp2
	call getarg(7,temparg); read(temparg,"(I4)") species2(1)
	call getarg(8,temparg); read(temparg,"(I4)") species2(2)
	call getarg(9,temparg); read(temparg,"(I4)") species2(3)
	call getarg(10,temparg); read(temparg,"(F12.4)") framestep
	call getarg(11,temparg); read(temparg,"(I4)") sp1
	nsites = (nargs - 11) / 2
	if (nsites.EQ.0) stop "No hydrogen bonding sites defined!"
	if (mod((nargs-11),2).EQ.1) stop "Incomplete site spec given."
	allocate(sites(nsites,2))
	do n=1,nsites
	  call getarg((2*n)+10,temparg); read(temparg,"(I4)") sites(n,1)
	  call getarg((2*n)+11,temparg); read(temparg,"(I4)") sites(n,2)
	end do
	
	write(0,*) "O...site distance min/max :",minoo,maxoo
	write(0,*) "Minimum O-H...site angle  :",minang

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	allocate(life(s_nmols(sp2),nsites))
	allocate(hist(framestodo))

	!sites(:,:) = reshape( (/ 6,8,9,10,11,12,0,20,21,22,23,24 /) , (/ nsites,2 /) )
	life = 0
	hist = 0.0
	nhbonds = 0

	write(0,"(A,4I3)") "'Water' species, X, H, H :",sp2,species2
	write(0,*) "XH site specification (X, H):"
	do m=1,nsites
	  write(0,*) (sites(m,n),n=1,2)
	end do

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

	i = s_start(sp1)-1
	!do m1=1,s_nmols(sp1)

	j = s_start(sp2)-1
	do m2=1,s_nmols(sp2)
	  c(1)=xpos(j+species2(1))	! Species2 electronegative site (e.g. oxygen in water)
	  c(2)=ypos(j+species2(1))
	  c(3)=zpos(j+species2(1))
	  do m=1,nsites
	    a(1)=xpos(sites(m,1))	! Species1 electronegative site (e.g. O in glucose OH)
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
		life(m2,m) = life(m2,m) + 1
		nhbonds = nhbonds + 1
		do o=1,life(m2,m)
		  hist(o) = hist(o) + 1
		end do
	      else
		! Reset lifetime counter
		life(m2,m) = 0
	      end if
	    end if
	    ! Now for H-bonds where the glucose O is the acceptor (glucO...waterH-waterO)
	    if (species2(2).NE.0) then
	      write(0,*) "Not yet written!"
	      ! For multiple ways for sp2 to H-bond, need to consider individual types of interaction as well (donor/acceptor)
	      b(1)=xpos(j+species2(2))	! Water Hydrogen 1
	      b(2)=ypos(j+species2(2))
	      b(3)=zpos(j+species2(2))
	      ! gO-wO distance (calc again to be safe)
	      call pbc(a(1),a(2),a(3),c(1),c(2),c(3),temp(1),temp(2),temp(3))
	      roo=sqrt( (temp(1)-c(1))**2 + (temp(2)-c(2))**2 + (temp(3)-c(3))**2 )
	      ! gO-wH-wO angle
              angle = calcangle(a,b,c)
	      if ((roo.LT.maxoo).AND.(angle.GT.minang)) then
		nhbonds = nhbonds + 1
		life(n,m) = life(n,m) + 1
		do o=1,life(n,m)
		  hist(o) = hist(o) + 1
		end do
	      else
		! Reset lifetime counter
		life(n,m) = 0
	      end if
	      b(1)=xpos(j+species2(3))	! Water Hydrogen 2
	      b(2)=ypos(j+species2(3))
	      b(3)=zpos(j+species2(3))
	      ! O-O Dist is still the same, so just need angle.
              angle = calcangle(a,b,c)
	      if ((roo.LT.maxoo).AND.(angle.GT.minang)) then
		nhbonds = nhbonds + 1
		life(n,m) = life(n,m) + 1
		do o=1,life(n,m)
		  hist(o) = hist(o) + 1
		end do
	      else
		! Reset lifetime counter
		life(n,m) = 0
	      end if
	    end if
	  end do
	  j = j+s_natoms(sp2)
	end do
	      
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


	open(unit=9,file="lifehb.dat",form="formatted",status="replace")

	! Normalise and write out the histogram
	do n=1,nframes
	  hist(n) = hist(n) / real(nhbonds)
	  write(9,"(3F10.4)") n*tstep, hist(n)
	end do
	close(9)

	write(0,*) "Finished."
999	close(14)
	end program lifehb

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
