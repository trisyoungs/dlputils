!	** vac **
!	Calculate the velocity autocorrelation

	program vacf
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	integer :: i,n,m,nframes,success,nargs,baselen,s
	integer :: framestodo, count1
	integer :: iargc
	real*8, allocatable :: vx(:,:), vy(:,:), vz(:,:), vac(:,:)
	real*8 :: tx,ty,tz,sp,deltat

	nargs = iargc()
	if (nargs.NE.4) then
	  write(0,"(A)") "Usage : vacf <HISTORYfile> <OUTPUTfile> <delta t> <framestodo>"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(F10.6)") deltat
	call getarg(4,temp); read(temp,"(I6)") framestodo

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
	 ! Now, read in the history header so that we have cell()
	 if (readheader().EQ.-1) goto 799

	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....

	allocate (vac(0:framestodo,natms))
	allocate (vx(framestodo,natms))
	allocate (vy(framestodo,natms))
	allocate (vz(framestodo,natms))

	vac = 0.0
	vx = 0.0
	vy = 0.0
	vz = 0.0
	
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	! Store all velocities
	vx(nframes,:) = xvel(:)
	vy(nframes,:) = yvel(:)
	vz(nframes,:) = zvel(:)

	! Next frame
	if (nframes.ne.framestodo) goto 101

	! Calculate VACF
	write(0,*) "Calculating VACF..."

	do n=1,nframes

	  write(0,"(i6,a,i6)") n,"/",nframes

	  do m=n,nframes

	    do i=1,natms

	      ! Calculate VACF to time t0(n) -> t1(m)
	      ! Calculate scalar product
	      sp = vx(n,i)*vx(m,i) + vy(n,i)*vy(m,i) + vz(n,i)*vz(m,i)
	      vac(m-n,i) = vac(m-n,i) + sp

	    end do

	  end do 

	end do

	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended before framestodo was fulfilled..."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",framestodo,")"
	goto 801
800	write(0,*) "Framestodo was fulfilled."
801	write(0,*) ""

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="vacresults."
	  baselen=11
	else
	  basename=hisfile(1:baselen)
	endif

	! Open output files
	do s=1,nspecies
   
	  resfile=basename(1:baselen)//"vac"//CHAR(48+s)
	  open(unit=20+s,file=resfile,form="formatted",status="replace")

	end do

	! Average over 'origins'
	do n=0,nframes-1

	  do i=1,natms

	    vac(n,i) = vac(n,i) / real(nframes - n)

	  end do

	end do

	! Average atomic VACF over molecules

	do n=0,nframes-1

	  do s=1,nspecies

	    do m=2,s_nmols(s)

	      do i=1,s_natoms(s)

	        ! Sum VACFs into that for first molecule
	        vac(n,s_start(s)+i-1) = vac(n,s_start(s)+i-1) + vac(n,s_start(s)+(m-1)*s_natoms(s)+i-1)

	      end do

	    end do

	    ! Average over number of molecules
	    do i=1,s_natoms(s)
	      vac(n,s_start(s)+i-1) = vac(n,s_start(s)+i-1) / s_nmols(s)
	    end do
	
	  end do

	  ! Write the data...
	  do s=1,nspecies
	    write(20+s,"(20F10.4)") n*deltat,(vac(n,s_start(s)+i-1),i=1,s_natoms(s))
	  end do

	end do

	do n=1,nspecies
	  close(20+n)
	end do

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program vacf

