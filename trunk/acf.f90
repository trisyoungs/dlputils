!	** acf **
!	Calculate the specified autocorrelation tensor

	program acf
	use dlprw; use utility
	implicit none
	include "mpif.h"

	integer, parameter :: VELOCITY=1, DIPOLE=2, FROMFILE=3
	logical :: MASTER, SLAVE
	character*80 :: hisfile1, dlpoutfile,basename,resfile,headerfile
	character*20 :: temp
	character*4 :: namepart
	integer :: i,j,n,m,nframes,success,nargs,baselen,sp,pos, t0, tn
	integer :: length, acftype = 0, framestodo = -1
	integer :: iargc, idx
	integer :: t_start, t_finish, err_mpi, id_mpi, nproc_mpi, qmax
	real*8, allocatable :: acf_total(:), acf_intra(:), accum(:), acfpart_total(:,:), acfpart_intra(:,:), qx(:,:), qy(:,:), qz(:,:)
	real*8, allocatable :: temp_total(:), temp_intra(:), temppart_total(:,:), temppart_intra(:,:), tempaccum(:)
	real*8 :: deltat, totmass, dist, vec(3)
	real*8 :: xx, yy, zz, xy, xz, yz
	logical :: altheader = .false.

	nargs = iargc()
	if (nargs.ne.7) then
	  write(0,"(a)") "Usage : acf <HISfile|quantityfile> <DLP OUTPUTfile> <headerfile [0 for none]> <functype=velocity,dipole,file> <delta t> <length> <nframes>"
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
	if (temp.eq."file") acftype = FROMFILE
	if (acftype.eq.0) stop "Invalid correlation function requested - options are 'velocity' or 'dipole'."
	call getarg(5,temp); read(temp,"(F10.6)") deltat
	call getarg(6,temp); read(temp,"(I6)") length
	call getarg(7,temp); read(temp,"(I6)") framestodo
	
	if (acftype.eq.VELOCITY) write(0,*) "Calculating velocity autocorrelation function."
	if (acftype.eq.DIPOLE) write(0,*) "Calculating molecular dipole autocorrelation function."
	if (acftype.eq.FROMFILE) write(0,*) "Calculating autocorrelation function from quantities in file."

	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
	qmax = sum(s_nmols)

	! Ascertain length of basename....
	baselen=-1
	do i=80,1,-1
	  if (hisfile1(i:i).eq.".") then
	    baselen=i
	    goto 50
	  endif
	end do
50     if (baselen.EQ.-1) then
	  basename="acfresults."
	  baselen=11
	else
	  basename=hisfile1(1:baselen)
	endif

	! Get part of output filename
	if (acftype.eq.VELOCITY) namepart = "vacf"
	if (acftype.eq.DIPOLE) namepart = "dacf"
	if (acftype.eq.FROMFILE) namepart = "facf"

	allocate (acf_total(0:length-1))
	allocate (acf_intra(0:length-1))
	allocate (acfpart_total(6,0:length-1))
	allocate (acfpart_intra(6,0:length-1))
	allocate (accum(0:length))
	allocate (qx(length,qmax))
	allocate (qy(length,qmax))
	allocate (qz(length,qmax))

	acf_total = 0.0
	acf_intra = 0.0
	acfpart_intra = 0.0
	acfpart_total = 0.0
	accum = 0.0
	qx = 0.0
	qy = 0.0
	qz = 0.0

	! Check history file (if not FROMFILE) or open binary quantity file
	if (acftype.NE.FROMFILE) then
60	  write(0,*) "Reading history file...."
 	  if (altheader) then
	    call openhis(headerfile,10)
	    if (readheader().EQ.-1) goto 799
	    close(dlpun_his)
	    call openhis(hisfile1,10)
	  else 
	    call openhis(hisfile1,10)
	   if (readheader().EQ.-1) goto 799
	  end if
	else
	  open(unit=11,file=hisfile1,form='unformatted',status='old')
	end if

	nframes=0
	pos=0

	! Initialise MPI and determine slave data
	call MPI_Init(err_mpi)
	call MPI_Comm_rank(MPI_COMM_WORLD,id_mpi,err_mpi)
	call MPI_Comm_size(MPI_COMM_WORLD,nproc_mpi,err_mpi)

	if ((nproc_mpi.le.1).or.(id_mpi.eq.0)) then
	  MASTER = .true.; SLAVE = .false.
	else
	  MASTER = .false.; SLAVE = .true.
	end if

	! Allocate extra arrays for MASTER
	if (MASTER) then
	  allocate (temp_total(0:length-1))
	  allocate (temp_intra(0:length-1))
	  allocate (temppart_total(6,0:length-1))
	  allocate (temppart_intra(6,0:length-1))
	  allocate (tempaccum(0:length))
	end if

	! Calculate atom range for node
	t_start = id_mpi * (length / nproc_mpi)
	t_finish = t_start + (length / nproc_mpi) - 1
	if (id_mpi+1.eq.nproc_mpi) t_finish = t_finish + mod(length,length/nproc_mpi)
	write(13+id_mpi,"(a,i2,a,i7,a)") "Process ",id_mpi," thinks there are ",nproc_mpi," processes"
	write(13+id_mpi,"(a,i2,a,i7,a)") "Process ",id_mpi," thinks there are ",qmax," molecules in total"
	write(13+id_mpi,"(a,i2,a,i7,a)") "Process ",id_mpi," thinks there are ",length," points in the correlation array"
	write(13+id_mpi,"(a,i2,a,i7,a,i7,a)") "Process ",id_mpi," calculating intervals ",t_start," to ",t_finish,"..."


	! Calculate array position
101	pos=mod(nframes,length)+1
	nframes = nframes+1
	if (acftype.eq.FROMFILE) then
	  if (MASTER) then
	    read(11,end=102,err=102) qx(pos,:)
	    read(11,end=102,err=102) qy(pos,:)
	    read(11,end=102,err=102) qz(pos,:)
	    goto 105
102	    call MPI_BCast(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    goto 799  ! End of file encountered, or file error....
	    ! Send 'continue' flag
105	    call MPI_BCast(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	  else
	    ! Slaves must wait for the 'continue' flag
	    call MPI_BCast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_mpi)
	    if (i.eq.0) goto 799
	  end if
	  call MPI_BCast(qx(pos,:),qmax,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	  call MPI_BCast(qy(pos,:),qmax,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
	  call MPI_BCast(qz(pos,:),qmax,MPI_REAL8,0,MPI_COMM_WORLD,err_mpi)
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

	if (MASTER.and.mod(nframes,100).EQ.0) write(0,"(i)") nframes

	! Calculate quantities for correlation function (only if acftype != FROMFILE)
	if (acftype.ne.FROMFILE) then
	  qx(pos,:) = 0.0
	  qy(pos,:) = 0.0
	  qz(pos,:) = 0.0
	  idx = 1
	  do sp=1,nspecies
	    if (acftype.eq.VELOCITY) then
	     call calc_com
	     i = s_start(sp)
	     do m=1,s_nmols(sp)
	       totmass = 0.0
	       do j=i,i+s_natoms(sp)-1
	         totmass = totmass + mass(j)
	         qx(pos,idx) = qx(pos,idx) + mass(j)*xvel(j)
	         qy(pos,idx) = qy(pos,idx) + mass(j)*yvel(j)
	         qz(pos,idx) = qz(pos,idx) + mass(j)*zvel(j)
	       end do
	       !if (m.eq.1) write(0,*) "Frame ",pos, "comvx(mol1) = ",comvx(pos,1)
	       qx(pos,idx) = qx(pos,idx) / totmass
  	       qy(pos,idx) = qy(pos,idx) / totmass
	       qz(pos,idx) = qz(pos,idx) / totmass
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
	         qx(pos,idx) = qx(pos,idx) + charge(j)*vec(1)
	         qy(pos,idx) = qy(pos,idx) + charge(j)*vec(2)
	         qz(pos,idx) = qz(pos,idx) + charge(j)*vec(3)
	       end do
	       ! Convert dipole from q.angstrom to Debye (C.m)
	       qx(pos,idx) = qx(pos,idx) * 4.80321
	       qy(pos,idx) = qy(pos,idx) * 4.80321
	       qz(pos,idx) = qz(pos,idx) * 4.80321
	       idx = idx + 1
	       i = i + s_natoms(sp)
	     end do
	    endif
	  end do
	end if

	! Accumulate ACF.

	if (nframes.ge.length) then
	  ! write(0,*) "Accumulating"

	  ! Variables:
	  ! 'n' will represent the frame distance
	  ! 'pos' contains the array position of the last added point
	  ! 't0' is the t=0 point, with which all scalars are calcd
	  ! 'tn' is the nth point in time

	  t0 = pos+1
	  if (t0.gt.length) t0=t0-length

	  do n=t_start, t_finish

	    tn = t0+n
	    if (tn.gt.length) tn=tn-length

	    ! Calculate total scalar products

	    xx = sum(qx(t0,:))*sum(qx(tn,:))
	    yy = sum(qy(t0,:))*sum(qy(tn,:))
	    zz = sum(qz(t0,:))*sum(qz(tn,:))
	    xy = sum(qx(t0,:))*sum(qy(tn,:)) + sum(qy(t0,:))*sum(qx(tn,:))
	    xz = sum(qx(t0,:))*sum(qz(tn,:)) + sum(qz(t0,:))*sum(qx(tn,:))
	    yz = sum(qy(t0,:))*sum(qz(tn,:)) + sum(qz(t0,:))*sum(qy(tn,:))

	    accum(n) = accum(n) + 1.0

	    acf_total(n) = acf_total(n) + xx + yy + zz
	    ! Accumulate partial ACFs
	    acfpart_total(1,n) = acfpart_total(1,n) + xx
	    acfpart_total(2,n) = acfpart_total(2,n) + yy
	    acfpart_total(3,n) = acfpart_total(3,n) + zz
	    acfpart_total(4,n) = acfpart_total(4,n) + xy
	    acfpart_total(5,n) = acfpart_total(5,n) + xz
	    acfpart_total(6,n) = acfpart_total(6,n) + yz

	    ! Calculate molecular products

	    xx = sum(qx(t0,:)*qx(tn,:))
	    yy = sum(qy(t0,:)*qy(tn,:))
	    zz = sum(qz(t0,:)*qz(tn,:))
	    xy = sum(qx(t0,:)*qy(tn,:)) + sum(qy(t0,:)*qx(tn,:))
	    xz = sum(qx(t0,:)*qz(tn,:)) + sum(qz(t0,:)*qx(tn,:))
	    yz = sum(qy(t0,:)*qz(tn,:)) + sum(qz(t0,:)*qy(tn,:))

	    acf_intra(n) = acf_intra(n) + xx + yy + zz
	    ! Accumulate partial VACFs
	    acfpart_intra(1,n) = acfpart_intra(1,n) + xx
	    acfpart_intra(2,n) = acfpart_intra(2,n) + yy
	    acfpart_intra(3,n) = acfpart_intra(3,n) + zz
	    acfpart_intra(4,n) = acfpart_intra(4,n) + xy
	    acfpart_intra(5,n) = acfpart_intra(5,n) + xz
	    acfpart_intra(6,n) = acfpart_intra(6,n) + yz

	  end do

	  ! Write intermediate results file?
	  if (mod(nframes-length,100).EQ.0) then

	    ! Reduce acf data into temporary arrays on master
	    call MPI_Reduce(acf_total,temp_total, length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
	    call MPI_Reduce(acf_intra,temp_intra, length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
	    call MPI_Reduce(acfpart_total,temppart_total, 6*length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
	    call MPI_Reduce(acfpart_intra,temppart_intra, 6*length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
	    call MPI_Reduce(accum,tempaccum, length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)

	    if (MASTER) then
	
	      open(unit=20,file=basename(1:baselen)//namepart//CHAR(48+sp)//".total_temp",form="formatted",status="replace")
	      write(20,'("# SP,INTVL=",3i5)') sp
	      write(20,"(10a14)") "#t","total","xx","yy","zz","xy","xz","yz","acc"
	      open(unit=21,file=basename(1:baselen)//namepart//CHAR(48+sp)//".intra_temp",form="formatted",status="replace")
	      write(21,'("# SP,INTVL=",3i5)') sp
	      write(21,"(10a14)") "#t","total","xx","yy","zz","xy","xz","yz","acc"
	      open(unit=22,file=basename(1:baselen)//namepart//CHAR(48+sp)//".inter_temp",form="formatted",status="replace")
	      write(22,'("# SP,INTVL=",3i5)') sp
	      write(22,"(10a14)") "#t","total","xx","yy","zz","xy","xz","yz","acc"
	      do n=0,length-1
	        write(20,"(9f14.6,e14.6)") n*deltat, temp_total(n)/tempaccum(n), (temppart_total(m,n)/tempaccum(n),m=1,6), tempaccum(n)
	        write(21,"(9f14.6,e14.6)") n*deltat, temp_intra(n)/tempaccum(n), (temppart_intra(m,n)/tempaccum(n),m=1,6), tempaccum(n)
	        write(22,"(9f14.6,e14.6)") n*deltat, (temp_total(n)-temp_intra(n))/tempaccum(n), ((temppart_total(m,n)-temppart_intra(m,n))/tempaccum(n),m=1,6), tempaccum(n)
	      end do
	      close(20)
	      close(21)
	      close(22)
	    else
	    end if

	  end if

	endif

	if (nframes.eq.framestodo) goto 801
	goto 101

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",nframes,")"
	goto 801
801	write(0,*) ""

	! Gather acf data into temporary arrays on master
	call MPI_Reduce(acf_total,temp_total, length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
	call MPI_Reduce(acf_intra,temp_intra, length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
	call MPI_Reduce(acfpart_total,temppart_total, 6*length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
	call MPI_Reduce(acfpart_intra,temppart_intra, 6*length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
	call MPI_Reduce(accum,tempaccum, length, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,err_mpi)
	call MPI_Finalize(err_mpi)

	! Write final functions
	open(unit=20,file=basename(1:baselen)//namepart//CHAR(48+sp)//".total",form="formatted",status="replace")
	write(20,'("# SP,INTVL=",3i5)') sp
	write(20,"(10a14)") "#t","total","xx","yy","zz","xy","xz","yz","acc"
	open(unit=21,file=basename(1:baselen)//namepart//CHAR(48+sp)//".intra",form="formatted",status="replace")
	write(21,'("# SP,INTVL=",3i5)') sp
	write(21,"(10a14)") "#t","total","xx","yy","zz","xy","xz","yz","acc"
	open(unit=22,file=basename(1:baselen)//namepart//CHAR(48+sp)//".inter",form="formatted",status="replace")
	write(22,'("# SP,INTVL=",3i5)') sp
	write(22,"(10a14)") "#t","total","xx","yy","zz","xy","xz","yz","acc"
	do n=0,length-1
	  write(20,"(9f14.6,e14.6)") n*deltat, temp_total(n)/tempaccum(n), (temppart_total(m,n)/tempaccum(n),m=1,6), tempaccum(n)
	  write(21,"(9f14.6,e14.6)") n*deltat, temp_intra(n)/tempaccum(n), (temppart_intra(m,n)/tempaccum(n),m=1,6), tempaccum(n)
	  write(22,"(9f14.6,e14.6)") n*deltat, (temp_total(n)-temp_intra(n))/tempaccum(n), ((temppart_total(m,n)-temppart_intra(m,n))/tempaccum(n),m=1,6), tempaccum(n)
	end do
	close(20)
	close(21)
	close(22)
      
	write(0,*) "Finished."
999	close(10)
	close(13)
	end program acf
