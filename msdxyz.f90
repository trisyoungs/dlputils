!	** msdxyz **
!	Calculate the xyz components of the mean square displacement of molecules (COM)

	program msdxyz
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,m,o,s1,m1,nframes,success,nargs,sp1
	integer :: framestodo, count1, totmols
	real*8, allocatable :: lastx(:), lasty(:), lastz(:)
	real*8, allocatable :: spdisp_x(:,:),spdisp_y(:,:),spdisp_z(:,:)	! Cartesian molecule displacements of MSD
	real*8, allocatable :: msd_x(:,:),msd_y(:,:),msd_z(:,:)   		! Cartesian MSD per molecule
	real*8, allocatable :: avg_msd_x(:,:),avg_msd_y(:,:),avg_msd_z(:,:)   	! Final Cartesian MSD averaged over molecules
	real*8, allocatable :: dx(:), dy(:), dz(:)
	integer, allocatable :: orn(:,:)
	integer :: iargc
	real*8 :: tx,ty,tz,rij2,deltat

	nargs = iargc()
	if (nargs.LT.4) then
	  write(0,"(A)") "Usage : msdxyz <HISTORYfile> <OUTPUTfile> <delta t> <framestodo> [-axis sp x x y y]"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(F10.6)") deltat
	call getarg(4,temp); read(temp,"(I6)") framestodo

	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
	allocate(aa(nspecies,4))

	n = 4
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-axis")
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp1
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp1,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp1,2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp1,3)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(sp1,4)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Molecular axis for species ",sp1," calculated from: X=",aa(sp1,1),"->", &
		 & aa(sp1,2),", Y=0.5X->0.5(r(",aa(sp1,3),")->r(",aa(sp1,4),"))"
	      do m=1,3
		if ((axesAatoms(sp1,m).lt.1).or.(axesAatoms(sp1,m).gt.s_natoms(sp1))) stop "Atom id out of range for axes on this species!"
	      end do
	      axesAdefined(sp1) = .true.
	  end select
	end do

	! Open and check the files...
	call openhis(hisfile,10)
	! Now, read in the history header so that we have cell()
	if (readheader().EQ.-1) goto 799

	totmols = 0
	do n=1,nspecies
	  totmols = totmols + s_nmols(n)
	end do

	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....

	allocate (orn(framestodo,totmols))
	allocate (spdisp_x(framestodo,totmols))
	allocate (spdisp_y(framestodo,totmols))
	allocate (spdisp_z(framestodo,totmols))
	allocate (msd_x(framestodo,totmols))
	allocate (msd_y(framestodo,totmols))
	allocate (msd_z(framestodo,totmols))
	allocate (avg_msd_x(framestodo,nspecies))
	allocate (avg_msd_y(framestodo,nspecies))
	allocate (avg_msd_z(framestodo,nspecies))
	allocate (lastx(totmols))
	allocate (lasty(totmols))
	allocate (lastz(totmols))
	allocate (dx(totmols))
	allocate (dy(totmols))
	allocate (dz(totmols))

	orn = 0
	msd_x = 0.0
	msd_y = 0.0
	msd_z = 0.0
	avg_msd_x = 0.0
	avg_msd_y = 0.0
	avg_msd_z = 0.0
	spdisp_x = 0.0
	spdisp_y = 0.0
	spdisp_z = 0.0
	
	call alloc_axis
	! Store the initial positions
	count1 = 0
	do s1=1,nspecies
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1

	    lastx(count1) = axisox(s1,m1)
	    lasty(count1) = axisoy(s1,m1)
	    lastz(count1) = axisoz(s1,m1)

	  end do !m1
	end do !s1
	call genaxis

100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	call calc_com

	! Calculate XYZ displacements of species between this frame and the last
	count1 = 0
	do s1=1,nspecies
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1

	    call pbc(comx(s1,m1),comy(s1,m1),comz(s1,m1),lastx(count1),lasty(count1),lastz(count1),tx,ty,tz)
	    tx = lastx(count1) - tx
	    ty = lasty(count1) - ty
	    tz = lastz(count1) - tz

	    ! Rotate the displacement dx,dy,dx into the local frame of the molecule's previous position
	    dx(count1) = tx*axisx(sp1,m1,1) + ty*axisx(sp1,m1,2) + tz*axisx(sp1,m1,3)
	    dy(count1) = tx*axisy(sp1,m1,1) + ty*axisy(sp1,m1,2) + tz*axisy(sp1,m1,3)
	    dz(count1) = tx*axisz(sp1,m1,1) + ty*axisz(sp1,m1,2) + tz*axisz(sp1,m1,3)

	  end do !m1
	end do !s1

	! Store current positions for use in next iteration
	count1 = 0
	do s1=1,nspecies
	  do m1=1,s_nmols(s1)
	    count1 = count1 + 1

	    lastx(count1) = comx(s1,m1)
	    lasty(count1) = comy(s1,m1)
	    lastz(count1) = comz(s1,m1)

	  end do !m1
	end do !s1


	! Update averages over origins
	do n=1,nframes
	  count1 = 0
	  do s1=1,nspecies
	    do m1=1,s_nmols(s1)
		count1 = count1 + 1

		! Origin COM positions given by rxyz(m,count1)
		! 'To' COM position given by ryxz(nframes,count1)

		! Accumulate the incremental cartesian displacements of molecules(in the molecular axes)
		! Array is 'backwards' - i.e. low 'n' = largest accumulation (over most frames)
		spdisp_x(n,count1) = spdisp_x(n,count1) + dx(count1)
		spdisp_y(n,count1) = spdisp_y(n,count1) + dy(count1)
		spdisp_z(n,count1) = spdisp_z(n,count1) + dz(count1)

		! Increment MSD arrays
		msd_x(nframes-(n-1),count1) = msd_x(nframes-(n-1),count1) + spdisp_x(n,count1)**2
		msd_y(nframes-(n-1),count1) = msd_y(nframes-(n-1),count1) + spdisp_y(n,count1)**2
		msd_z(nframes-(n-1),count1) = msd_z(nframes-(n-1),count1) + spdisp_z(n,count1)**2
		orn(nframes-(n-1),count1) = orn(nframes-(n-1),count1) + 1

	    end do !m1
	  end do !s1
	end do !frames

	call genaxis

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

	open(unit=9,file="msd.dat",form="formatted",status="replace")
	do n=1,nframes
	  count1 = 0
	  do s1=1,nspecies
	    do m1=1,s_nmols(s1)
		count1 = count1 + 1

		! Sum over number of molecules
		avg_msd_x(n,s1) = avg_msd_x(n,s1) + msd_x(n,count1)
		avg_msd_y(n,s1) = avg_msd_y(n,s1) + msd_y(n,count1)
		avg_msd_z(n,s1) = avg_msd_z(n,s1) + msd_z(n,count1)

	    end do

	    ! Average over number of molecules and origins
	    avg_msd_x(n,s1) = avg_msd_x(n,s1) / real(orn(n,s1)) / s_nmols(s1)
	    avg_msd_y(n,s1) = avg_msd_y(n,s1) / real(orn(n,s1)) / s_nmols(s1)
	    avg_msd_z(n,s1) = avg_msd_z(n,s1) / real(orn(n,s1)) / s_nmols(s1)
	
	  end do

	  ! Write the data...
	  write(9,"(10F10.4)") n*deltat,((avg_msd_x(n,m),avg_msd_y(n,m),avg_msd_z(n,m),m=1,nspecies),o=1,3)
	end do
	close(9)

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program msdxyz

