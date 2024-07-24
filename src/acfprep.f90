!	** acfprep **
!	Calculate the specified quantity ready for use with ACF

	program acfprep
	use dlprw; use utility
	implicit none

	integer, parameter :: VELOCITY=1, DIPOLE=2
	character*80 :: hisfile1, dlpoutfile,basename,resfile,headerfile
	character*20 :: temp
	character*4 :: namepart
	integer :: i,j,n,m,nframes,success,nargs,baselen,sp,idx
	integer :: length, acftype = 0, framestodo = -1
	integer :: iargc, qmax
	real*8, allocatable :: qx(:), qy(:), qz(:)
	real*8 :: totmass, dist, vec(3)
	real*8 :: xx, yy, zz, xy, xz, yz
	logical :: altheader = .false.

	nargs = iargc()
	if (nargs.ne.5) then
	  write(0,"(a)") "Usage : acfprep <HISTORYfile> <OUTPUTfile> <headerfile [0 for none]> <functype=velocity,dipole> <outputfile>"
	  stop
	end if
	call getarg(1,hisfile1)
	call getarg(2,dlpoutfile)
	call getarg(3,headerfile)
	if (headerfile.ne."0") then
	  altheader = .true.
	  write(0,*) "Alternative header file supplied."
	end if
	call getarg(4,temp)
	if (temp.eq."velocity") acftype = VELOCITY
	if (temp.eq."dipole") acftype = DIPOLE
	if (acftype.eq.0) stop "Invalid quantity requested - options are 'velocity' or 'dipole'."
	call getarg(5,resfile)
	
	if (acftype.eq.VELOCITY) write(0,*) "Calculating COM velocities."
	if (acftype.eq.DIPOLE) write(0,*) "Calculating molecular dipoles."

	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
	qmax = sum(s_nmols)
	write(0,*) "Total number of molecules is ", qmax

	allocate (qx(qmax))
	allocate (qy(qmax))
	allocate (qz(qmax))

	qx = 0.0
	qy = 0.0
	qz = 0.0

	! Check history file
60	write(0,*) "Reading history file...."
 	if (altheader) then
	  call openhis(headerfile,10)
	  if (readheader().EQ.-1) goto 799
	 close(dlpun_his)
	  call openhis(hisfile1,10)
	else 
	  call openhis(hisfile1,10)
	 if (readheader().EQ.-1) goto 799
	end if

	! Open output file (binary)
	open(unit=15,file=resfile,form='unformatted',status='new')

	nframes=0

101	success=readframe()
	if (success.ne.0) goto 799  ! End of file encountered, or file error....

	nframes = nframes+1
	if (mod(nframes,100).EQ.0) write(0,"(i10)") nframes

	! Calculate quantities for correlation function, writing after each species
	idx = 1
	qx = 0.0
	qy = 0.0
	qz = 0.0
	do sp=1,nspecies
	  if (acftype.eq.VELOCITY) then
	    call calc_com
	    i = s_start(sp)
	    do m=1,s_nmols(sp)
	      totmass = 0.0
	      do j=i,i+s_natoms(sp)-1
	        totmass = totmass + mass(j)
	        qx(idx) = qx(idx) + mass(j)*xvel(j)
	        qy(idx) = qy(idx) + mass(j)*yvel(j)
	        qz(idx) = qz(idx) + mass(j)*zvel(j)
	      end do
	      qx(idx) = qx(idx) / totmass
	      qy(idx) = qy(idx) / totmass
	      qz(idx) = qz(idx) / totmass
	      idx = idx + 1
	      i = i + s_natoms(sp)
	    end do
	  else if (acftype.eq.DIPOLE) then
	    ! Take all atom positions relative to first atom in molecule
	    i = s_start(sp)
	    do m=1,s_nmols(sp)
	      do j=i,i+s_natoms(sp)-1
	        ! Get mim position of this atom with first
	        call pbc(xpos(j),ypos(j),zpos(j),xpos(i),ypos(i),zpos(i),vec(1),vec(2),vec(3))
	        qx(idx) = qx(idx) + charge(j)*vec(1)
	        qy(idx) = qy(idx) + charge(j)*vec(2)
	        qz(idx) = qz(idx) + charge(j)*vec(3)
	      end do
	      idx = idx + 1
	      i = i + s_natoms(sp)
	    end do
	  end if
	end do

	! Convert dipole from q.angstrom to Debye (C.m)
	if (acftype.eq.DIPOLE) then
	  qx = qx * 4.80321d0
	  qy = qy * 4.80321d0
	  qz = qz * 4.80321d0
	end if

	! Write data for all molecules
	write(15) qx
	write(15) qy
	write(15) qz

	if (nframes.eq.framestodo) goto 801
	goto 101

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended"
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",nframes,")"
	goto 801
801	write(0,*) ""

	close(15)

999	write(0,*) "Finished."
	end program acfprep
