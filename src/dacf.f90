!	** dacf **
!	Calculate the system dipole ACF

	program dacf
	use dlprw; use utility
	implicit none
	include "mpif.h"

	integer, parameter :: DIPOLE=2, FROMFILE=3
	logical :: MASTER, SLAVE
	character*80 :: hisfile1, dlpoutfile,basename,resfile,headerfile
	character*20 :: temp
	character*4 :: namepart
	integer :: i,j,n,m,nframes,success,nargs,baselen,sp,pos, t0, tn
	integer :: length, acftype = 0, framestodo = -1, framestodiscard = -1
	integer :: iargc, idx
	integer :: t_first, t_last, err_mpi, id_mpi, nproc_mpi, qmax
	real*8, allocatable :: acf_total(:), accum(:), qx(:,:), qy(:,:), qz(:,:), qtot(:,:)
	real*8, allocatable :: qxcurrent(:), qycurrent(:), qzcurrent(:)
	real*8, allocatable :: temp_total(:), tempaccum(:)
	real*8 :: deltat, dist, vec(3), qtotcurrent(3)
	real*8 :: xx, yy, zz
	logical :: altheader = .false.

	nargs = iargc()
	if (nargs.lt.7) then
	  write(0,"(a)") "Usage : acf <HISfile|quantityfile> <OUTPUTfile> <headerfile [0 for none]> <functype=velocity,dipole,file> <delta t> <length> <framestodo> [skip]"
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
	if (temp.eq."dipole") acftype = DIPOLE
	if (temp.eq."file") acftype = FROMFILE
	if (acftype.eq.0) stop "Invalid correlation function requested - options are 'velocity' or 'dipole'."
	call getarg(5,temp); read(temp,"(F10.6)") deltat
	call getarg(6,temp); read(temp,"(I10)") length
	call getarg(7,temp); read(temp,"(I10)") framestodo
	if (nargs.eq.8) then
	  call getarg(8,temp); read(temp,"(I10)") framestodiscard
	end if

	! Initialise MPI and determine slave data
	call MPI_Init(err_mpi)
	call MPI_Comm_rank(MPI_COMM_WORLD,id_mpi,err_mpi)
	call MPI_Comm_size(MPI_COMM_WORLD,nproc_mpi,err_mpi)

	if ((nproc_mpi.le.1).or.(id_mpi.eq.0)) then
	  MASTER = .true.; SLAVE = .false.
	else
	  MASTER = .false.; SLAVE = .true.
	end if

	if (MASTER) then
	  if (acftype.eq.DIPOLE) write(0,*) "Calculating molecular dipole autocorrelation function."
	  if (acftype.eq.FROMFILE) write(0,*) "Calculating autocorrelation function from quantities in file."
	  if (outinfo(dlpoutfile,1).EQ.-1) then
	    write(0,*) "Problem with OUTPUT file."
	    ! Send 'failed' flag
	    call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    goto 801
	  end if
	  ! Send 'continue' flag and pass necessary data to slaves
	  call MPI_BCast(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  call MPI_BCast(nspecies,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  call MPI_BCast(s_nmols,nspecies,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  call MPI_BCast(s_natoms,nspecies,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  call MPI_BCast(s_start,nspecies,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  ! Ascertain length of basename....
	  baselen=-1
	  do i=80,1,-1
	    if (hisfile1(i:i).eq.".") then
	      baselen=i
	      goto 50
	    endif
	  end do
50	  if (baselen.EQ.-1) then
	    basename="acfresults."
	    baselen=11
	  else
	    basename=hisfile1(1:baselen)
	  endif

	  ! Get part of output filename
	  if (acftype.eq.DIPOLE) namepart = "dacf"
	  if (acftype.eq.FROMFILE) namepart = "facf"

	else
	  ! Slaves must wait for the 'continue' flag
	  call MPI_BCast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  if (i.eq.0) goto 801
	  call MPI_BCast(nspecies,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  allocate(s_nmols(nspecies),s_natoms(nspecies),s_start(nspecies))
	  call MPI_BCast(s_nmols,nspecies,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  call MPI_BCast(s_natoms,nspecies,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  call MPI_BCast(s_start,nspecies,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	end if
	qmax = sum(s_nmols)
	  

	allocate (acf_total(0:length-1))
	allocate (accum(0:length))

	! Calculate time range for node and assign arrays
	t_first = id_mpi * (length / nproc_mpi) + 1
	t_last = t_first + (length / nproc_mpi) - 1
	if (id_mpi+1.eq.nproc_mpi) t_last = t_last + mod(length,length/nproc_mpi)
	write(13+id_mpi,"(a,i2,a,i7,a)") "Process ",id_mpi," thinks there are ",nproc_mpi," processes"
	write(13+id_mpi,"(a,i2,a,i7,a)") "Process ",id_mpi," thinks there are ",qmax," molecules in total"
	write(13+id_mpi,"(a,i2,a,i7,a)") "Process ",id_mpi," thinks there are ",length," points in the correlation array"
	write(13+id_mpi,"(a,i2,a,i7,a,i7,a)") "Process ",id_mpi," stores intervals ",t_first," to ",t_last,"..."

	allocate (qx(t_first:t_last,qmax))
	allocate (qy(t_first:t_last,qmax))
	allocate (qz(t_first:t_last,qmax))
	allocate (qtot(t_first:t_last,3))
	allocate (qxcurrent(qmax))
	allocate (qycurrent(qmax))
	allocate (qzcurrent(qmax))

	acf_total = 0.0
	accum = 0.0
	qx = 0.0
	qy = 0.0
	qz = 0.0
	qtot = 0.0
	qxcurrent = 0.0
	qycurrent = 0.0
	qzcurrent = 0.0
	qtotcurrent = 0.0

	if (MASTER) then
	  ! Check history file (if not FROMFILE) or open binary quantity file
	  if (acftype.NE.FROMFILE) then
60	    write(0,*) "Reading history file...."
 	    if (altheader) then
	      call openhis(headerfile,10)
	      if (readheader().EQ.-1) then
		! Send 'fail' flag
		call MPI_BCast(799,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
		goto 799
	      end if
	      close(dlpun_his)
	      call openhis(hisfile1,10)
	    else 
	      call openhis(hisfile1,10)
	     if (readheader().EQ.-1) then
		! Send 'fail' flag
		call MPI_BCast(799,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
		goto 799
	      end if
	    end if
	  else
	    open(unit=11,file=hisfile1,form='unformatted',status='old',err=65)
	    goto 70
	    ! Send 'fail' flag
65	    call MPI_BCast(798,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    write(0,*) "Couldn't open quantity file."
70	  end if
	  ! Send 'continue' flag
	  call MPI_BCast(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	else
	  ! Slaves must wait for the 'continue' flag
	  call MPI_BCast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  if (i.eq.799) goto 799
	  if (i.eq.798) goto 798
	end if

	write(0,*) "Process is ready", id_mpi
	nframes=0
	pos=0

	! Allocate extra arrays for MASTER
	if (MASTER) then
	  allocate (temp_total(0:length-1))
	  allocate (tempaccum(0:length))
	end if

	! Calculate array position
101	pos=mod(nframes,length)+1
	if (acftype.eq.FROMFILE) then
	  if (MASTER) then
	    read(11,end=102,err=102) qxcurrent
	    read(11,end=102,err=102) qycurrent
	    read(11,end=102,err=102) qzcurrent
	    qtotcurrent(1) = sum(qxcurrent)
	    qtotcurrent(2) = sum(qycurrent)
	    qtotcurrent(3) = sum(qzcurrent)
	!write(22,"(a,i3,a,3f10.4)") "n = ",pos, " q=",qtotcurrent(:)
	    goto 105
102	    call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    goto 798  ! End of file encountered, or file error....
	    ! Send 'continue' flag
105	    call MPI_BCast(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  else
	    ! Slaves must wait for the 'continue' flag
	    call MPI_BCast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    if (i.eq.0) goto 798
	  end if
	  call MPI_BCast(qtotcurrent,3,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	else
	  if (MASTER) then
	    success=readframe()
	    if (success.ne.0) then
	      call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	      goto 799  ! End of file encountered, or file error....
	    end if
	    ! Send 'continue' flag
	    call MPI_BCast(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  else
	    ! Slaves must wait for the 'continue' flag
	    call MPI_BCast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    if (i.eq.0) goto 799
	  end if
	  ! Broadcast frame data
	  call MPI_BCast(xpos,natms,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	  call MPI_BCast(ypos,natms,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	  call MPI_BCast(zpos,natms,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	  ! Send charge and mass data on frame 1
	  if (nframes.eq.1) then
	    call MPI_BCast(charge,natms,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	    call MPI_BCast(mass,natms,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	  end if
	end if

	if (framestodiscard.gt.0) then
	  if (MASTER.and.mod(framestodiscard,100).EQ.0) write(0,"('Skipping ',i6)") framestodiscard
	  framestodiscard = framestodiscard - 1
	  goto 101
	end if

	nframes = nframes+1
	if (MASTER.and.mod(nframes,100).EQ.0) write(0,"(i10)") nframes

	! Calculate quantities for correlation function (only if acftype != FROMFILE)
	if (acftype.ne.FROMFILE) then
	  qxcurrent = 0.0
	  qycurrent = 0.0
	  qzcurrent = 0.0
	  idx = 1
	  do sp=1,nspecies
	   ! Take all atom positions relative to first atom in molecule
	   i = s_start(sp)
	   do m=1,s_nmols(sp)
	     do j=i,i+s_natoms(sp)-1
	       ! Get mim position of this atom with first
	       call pbc(xpos(j),ypos(j),zpos(j),xpos(i),ypos(i),zpos(i),vec(1),vec(2),vec(3))
	       qxcurrent(idx) = qxcurrent(idx) + charge(j)*vec(1)
	       qycurrent(idx) = qycurrent(idx) + charge(j)*vec(2)
	       qzcurrent(idx) = qzcurrent(idx) + charge(j)*vec(3)
	     end do
	     idx = idx + 1
	     i = i + s_natoms(sp)
	   end do
	  end do
	  qtotcurrent(1) = sum(qxcurrent)
	  qtotcurrent(2) = sum(qycurrent)
	  qtotcurrent(3) = sum(qzcurrent)
	  ! Convert dipole from q.angstrom to Debye (C.m)
	  qtotcurrent = qtotcurrent * 4.80321
	end if

	! Store new quantity data, but only if falls within the frame range stored by this process
	if ((pos.ge.t_first).and.(pos.le.t_last)) then
	  qtot(pos,:) = qtotcurrent(:)
	end if

	! Now, determine new origin and who owns, and then send it out to all processes
	t0 = pos + 1
	if (t0.gt.length) t0=t0-length
	if ((t0.ge.t_first).and.(t0.le.t_last)) then
	  qtotcurrent(:) = qtot(t0,:)
	end if
	! Work out who has the data... (take care, since the last process may have slightly more data than the others)
	i = (t0-1) / (length/nproc_mpi)
	if (i.ge.nproc_mpi) i = nproc_mpi-1
	!write(0,*) "At t0=",t0, "Process ",i,"has the next qtotcurrent", length,nproc_mpi
	call MPI_BCast(qtotcurrent,3,MPI_REAL8,i,MPI_COMM_WORLD,err_mpi)

	! Accumulate ACF.
	if (nframes.ge.length) then

	  ! Variables:
	  ! 'n' will represent the correlation 'distance'
	  ! 't0' contains the current time origin that we're considering
	  ! 'tn' is the nth point in time available on this slave

	  do tn=t_first,t_last

	    ! Determine distance between this data 'tn' and the 'current' quantity data, wrapping value if necessary
	    n = tn-t0
	    if (n.lt.0) n=n+length

	    ! Calculate total system dipole autocorrelation function

	    acf_total(n) = acf_total(n) + dot_product(qtotcurrent(:),qtot(tn,:))
	!write(0,"(a,i3,a,3i3,a,2i3,a,3f10.4,a,3f10.4,a,f10.4)") "Process ",id_mpi," tn/tfirst/tlast=",tn,t_first,t_last," t0,n=",t0,n," qcrnt=",qtotcurrent(:)," qtot=",qtot(tn,:)," dp=",dot_product(qtotcurrent(:),qtot(tn,:))
	    accum(n) = accum(n) + 1.0

	  end do

	  framestodo = framestodo -1

	endif

	if (framestodo.eq.0) goto 801
	goto 101

798	if (MASTER) then
	  write(0,*) "Quantity file ended."
	  write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",nframes,")"
	end if
	goto 801
799	if (MASTER) then
	  write(0,*) "HISTORY file ended."
	  write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",nframes,")"
	end if
801	write(0,*) "Done", id_mpi

	! Gather acf data into temporary arrays on master
	call MPI_Reduce(acf_total,temp_total, length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
	call MPI_Reduce(accum,tempaccum, length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
 
	if (MASTER) then
	  ! Write final functions
	  open(unit=20,file=basename(1:baselen)//namepart//".total",form="formatted",status="replace")
	  !write(20,"(10a14)") "#t","total","xx","yy","zz","xy","xz","yz","acc"
	  do n=0,length-1
	    write(20,"(2f14.6,e14.6)") n*deltat, temp_total(n)/tempaccum(n), tempaccum(n)
	  end do
	  close(20)
	end if

	call MPI_Finalize(err_mpi)
      
	write(0,*) "Finished."
999	close(10)
	close(13)
	end program dacf
