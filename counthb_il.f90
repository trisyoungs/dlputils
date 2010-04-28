!	** counthb **
!	Counts the number of h-bonds between specified species per frame.

	program counthb_il
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,resfile
	character*20 :: temparg
	integer :: n,nframes,success,nargs,m,i,cat,an,j
	integer :: framestodo, nspec, t(9)
	integer :: anionsp = 2, numch = 9, cationsp = 1
	integer :: dgrid, agrid, nacceptor, ndonor, abin, dbin
	integer :: iargc
	real*8 :: a(3),b(3),c(3),temp(3),v1(3),v2(3),v1mag,v2mag,roo,angle,dp,maxoo,minang,totaverage,minoo
	integer, allocatable :: ilch(:,:)
	real*8, allocatable :: frameaverages(:), averages(:), results(:,:), totals(:)
	integer, parameter :: DONOR = 1, ACCEPTOR = 2
	real*8 :: ddelta = 0.01, adelta = 1.0

	minoo = 2.0

	nargs = iargc()
	if (nargs.LT.5) stop "Usage : counthb <DLP HISTORYfile> <DLP OUTPUTfile> <nframes> <maxoo> <minang>"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temparg); read(temparg,"(I5)") framestodo
	call getarg(4,temparg); read(temparg,"(F12.4)") maxoo
	call getarg(5,temparg); read(temparg,"(F12.4)") minang
	
	write(0,*) "C-Cl distance min/max :",minoo,maxoo
	write(0,*) "Minimum C-H-Cl angle  :",minang

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	allocate(results(framestodo,numch))
	allocate(ilch(numch,2))
	allocate(totals(framestodo))
	allocate(frameaverages(framestodo))
	allocate(averages(framestodo))
	ilch(:,:) = reshape( (/ 2,4,5,6,6,6,7,7,7,8,9,10,11,12,13,14,15,16 /) , (/ numch,2 /) )
	results(:,:) = reshape( (/ (0,n=1,numch*framestodo) /) , (/ framestodo,numch /) )
	totals = 0.0
	totaverage = 0.0
	frameaverages = 0.0
	averages = 0.0

	write(0,*) "Cation C-H specification:"
	do m=1,numch
	  write(0,*) (ilch(m,n),n=1,2)
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

	! Search for H-Bonds where cation C-H is the donor (C-H...O)
	i = s_start(cationsp)-1
	do cat=1,s_nmols(cationsp)		! Loop over cations
	  do n=1,numch				! Loop over C-H sites in cation
	    a(1)=xpos(i+ilch(n,1))	! Cation carbon
	    a(2)=ypos(i+ilch(n,1))
	    a(3)=zpos(i+ilch(n,1))
	    b(1)=xpos(i+ilch(n,2))	! Cation hydrogen
	    b(2)=ypos(i+ilch(n,2))
	    b(3)=zpos(i+ilch(n,2))
	    ! Get minimum image of atom a w.r.t. atom b
	    call pbc(a(1),a(2),a(3),b(1),b(2),b(3),temp(1),temp(2),temp(3))
	    a = temp
	    v1 = b - a
	    v1mag = sqrt( v1(1)**2 + v1(2)**2 + v1(3)**2 )
	    v1 = v1 / v1mag
	    ! Now loop over all chlorides
	    j = s_start(anionsp)
	    do an=1,s_nmols(anionsp)
	      c(1)=xpos(j)	! Chloride
	      c(2)=ypos(j)
	      c(3)=zpos(j)
	      ! Get minimum image of atom c w.r.t. atom b
	      call pbc(c(1),c(2),c(3),b(1),b(2),b(3),temp(1),temp(2),temp(3))
	      c = temp
	      ! Get C-Cl distance
	      roo = sqrt( (a(1)-c(1))**2 + (a(2)-c(2))**2 + (a(3)-c(3))**2 )
	      ! Calculate the angle
	      v2 = b - c
	      v2mag = sqrt( v2(1)**2 + v2(2)**2 + v2(3)**2 )
	      v2 = v2 / v2mag
              ! Calculate dot product and angle...
              dp = (v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3))
              angle = acos(dp)* 57.29577951d0
	      ! H-Bond?
	      if ((roo.LT.maxoo).AND.(angle.GT.minang)) then
		results(nframes,n) = results(nframes,n) + 1.0
	      end if
	      j = j + s_natoms(anionsp)
	    end do
	  end do
	  i = i+s_natoms(cationsp)
	end do
	      
	! Create the totals for this frame
	do n=1,numch
	  results(nframes,n) = results(nframes,n) / s_nmols(cationsp)
	  totals(nframes) = totals(nframes) + results(nframes,n)
	  averages(n) = averages(n) + results(nframes,n)
	end do
	totaverage = totaverage + totals(nframes)
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
	do n=1,numch
	  averages(n) = averages(n) / nframes
	end do
	
	open(unit=9,file="counthb.dat",form="formatted",status="replace")
	write(9,"(A)") "Frame   H2    H4   H5    H11   H12   H13   H14   H15   H16   Total"
	do n=1,framestodo
	  write(9,"(I5,2x,13(F5.3,1X),F7.4)") n,(results(n,m),m=1,numch),totals(n),frameaverages(n)
	end do
	write(9,"(' Avg ',2x,13(F5.3,1X),F7.4)") (averages(n),n=1,numch),totaverage
	close(9)

	write(0,*) "Finished."
999	close(14)
	end program counthb_il

