	! Compare molecules in crystal structure reference configuration with those in MD trajectory frames

	program cryscomp
	use parse; use dlprw; use utility
	implicit none
	character*80 :: hisfile, outfile, configfile, basename, altheaderfile, resfile
	character*20 :: temp
	integer :: nargs, n, sp, m, baselen, nframes, nframesused, success, startf = 1, endf = 0
	logical :: altheader = .false., npt = .false.
	real*8, allocatable :: refxyz(:,:,:), devxyz(:,:,:,:), spdevxyz(:,:,:), sprmsdxyz(:,:)
	real*8, allocatable :: refaxes(:,:,:), devangle(:,:,:,:), spdevangle(:,:,:), sprmsdangle(:,:)
	integer :: bin(4),nbins,nanglebins,tth,th,hun,units,ten
	real*8 :: tx, ty, tz, dist, dp(3), binwidth = 0.1, width = 5.0, anglebinwidth = 0.01
	integer :: iargc

        nargs = iargc()
        if (nargs.lt.2) then
	  write(0,"(a)") "Usage: cryscomp <HISTORYfile> <OUTfile> <CONFIGfile> ..."
	  write(0,"(a)") "        [-axis sp x1 x2 y1 y2]  Atoms to use for axis calculation in species sp"
	  write(*,"(a)") "        [-header file]          Use specified file to get header (for headerless HIS files)"
	  write(*,"(a)") "        [-start frameno]        Trajectory frame to start calculations (default = 1)"
	  write(*,"(a)") "        [-end frameno]          Trajectory frame to end calculations on (default = last)"
	  write(*,"(a)") "        [-binwidth r]           XYZ-displacement histogram bin width"
	  write(*,"(a)") "        [-width d]              XYZ-displacement histogram data range (-d:d) in Angstroms"
	  write(*,"(a)") "        [-npt d]                Assume NPT simulation and calculate relative displacements"
	  stop
	end if
        call getarg(1,hisfile)
        call getarg(2,outfile)
        call getarg(3,configfile)

	! Open files and get basename
	write(0,"(A,A)") "History file : ",hisfile
	write(0,"(A,A)") " Output file : ",outfile
	if (outinfo(outfile,1).eq.-1) stop "Can't read OUTPUT file."
	write(0,"(A,A)") " Config file : ",configfile

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  IF (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 50
	  endIF
	end do
50	if (baselen.eq.-1) then
	  basename="cryscomp."
	  baselen=9
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=15,file=basename(1:baselen)//"cc",form="formatted",status="replace")

	! Allocate arrays
	call alloc_axis()
	allocate(refxyz(nspecies,maxval(s_nmols),3))
	allocate(refaxes(nspecies,maxval(s_nmols),9))

	n = 3
        do
          n = n + 1; if (n.gt.nargs) exit
          call getarg(n,temp)
          select case (temp)
	    case ("-axis") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp,2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp,3)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp,4)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp," calculated from: X=",axesAatoms(sp,1),"->", &
	        & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
              do m=1,3
                if ((axesAatoms(sp,m).lt.1).or.(axesAatoms(sp,m).gt.s_natoms(sp))) stop "Atom id out of range for axes on this species!"
              end do
	      axesAdefined(sp) = .true.
	    case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
              altheader = .TRUE.
	    case ("-start")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") startf
	      write(0,"(A,I5)") "Starting frame = ",startf
	      write(15,"(A,I5)") "Starting frame = ",startf
	    case ("-end")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") endf
	      write(0,"(A,I5)") "End frame = ",endf
	      write(15,"(A,I5)") "End frame = ",endf
	    case ("-binwidth")
	      n = n + 1; call getarg(n,temp); read(temp,"(f8.5)") binwidth
	      write(0,"(A,f8.5)") "Histogram bin width = ",binwidth
	      write(15,"(A,f8.5)") "Histogram bin width = ",binwidth
	    case ("-width")
	      n = n + 1; call getarg(n,temp); read(temp,"(f8.5)") width
	      write(0,"(A,f8.5)") "Histogram data width (-d:d) = ",width
	      write(15,"(A,f8.5)") "Histogram data width (-d:d) = ",width
	    case default
	      write(0,*) "Unrecognised argument :",temp
	      stop
	  end select
	end do

	! Allocate more arrays
	nbins = width / binwidth
	nanglebins = 1.0 / anglebinwidth
	allocate(devxyz(nspecies,maxval(s_nmols),4,-nbins:nbins))
	allocate(spdevxyz(nspecies,4,-nbins:nbins))
	allocate(sprmsdxyz(nspecies,4))
	allocate(devangle(nspecies,maxval(s_nmols),4,-nanglebins:nanglebins))
	allocate(spdevangle(nspecies,4,-nanglebins:nanglebins))
	allocate(sprmsdangle(nspecies,4))

	! Print out a summary of the control variables to the output file.
	write(15,"(A,I3)") "Molecular species in file : ",nspecies
	do sp=1,nspecies
	  ! Check that molecules have had their axes defined
	  if (axesAdefined(sp)) then
	    write(15,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp," calculated from: X=",axesAatoms(sp,1),"->", &
	      & axesAatoms(sp,2),", Y=0.5X->0.5(r(",axesAatoms(sp,3),")->r(",axesAatoms(sp,4),"))"
	  else
	    if (s_natoms(sp).gt.2) stop "Axes must be defined for all molecular species where possible."
	  end if
	end do

	! Read the header of the history file...
        call openhis(hisfile,10)
        if (readheader().EQ.-1) then
          if (altheader) then
            write(0,*) "Restarted trajectory:"
            close(dlpun_his)
            call openhis(altheaderfile,10)
            if (readheader().EQ.-1) stop "Couldn't read header of alternative file."
            close(dlpun_his)
            call openhis(hisfile,10)
          else
            stop "Couldn't read header of history file."
          end if
        end if

	! Read in the configuration file and calculate our reference values
	n = natms
	call readconfig(configfile)
	! Check natoms is the same
	if (n.ne.natms) stop "Number of atoms in CONFIG file does not match number in HISTORY file!"
	! Transfer atom positions to work arrays and calculate axes
	do n=1,natms
	  xpos(n) = cfgxyz(n,1)
	  ypos(n) = cfgxyz(n,2)
	  zpos(n) = cfgxyz(n,3)
	end do
	call genaxes()
	call calc_com
	! Store axis systems and geometric centres
	do sp=1,nspecies
	  do m=1,s_nmols(sp)
	    refxyz(sp,m,1) = comx(sp,m)
	    refxyz(sp,m,2) = comy(sp,m)
	    refxyz(sp,m,3) = comz(sp,m)
	    refaxes(sp,m,:) = axesA(sp,m,:)
	  end do
	end do

	! Zero data before calculation
	devxyz = 0.0
	spdevxyz = 0.0
	sprmsdxyz = 0.0
	devangle = 0.0
	spdevangle = 0.0
	sprmsdangle = 0.0

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
	nframes = 0
	nframesused = 0
	! If we're using a mapfile, read it here
102	success=readframe()
	IF (success.eq.1) goto 120  ! End of file encountered....
	IF (success.lt.0) goto 119  ! File error....
	nframes=nframes+1
	if (MOD(nframes,100).eq.0) write(0,*) nframes
	if (nframes.LT.startf) goto 102

	nframesused = nframesused + 1

	! Generate all molecular axes and geometric centres
	call genaxes()
	call calc_com

	do sp=1,nspecies
	  do m=1,s_nmols(sp)

	    ! --------------------
	    ! Positional variation
	    ! --------------------

	    ! Calculate minimum image position of molecule m about 0,0,0
	    call pbc(comx(sp,m),comy(sp,m),comz(sp,m),refxyz(sp,m,1),refxyz(sp,m,2),refxyz(sp,m,3),tx,ty,tz)
	    tx = tx - refxyz(sp,m,1)
	    ty = ty - refxyz(sp,m,2)
	    tz = tz - refxyz(sp,m,3)
	    dist = sqrt(tx*tx + ty*ty + tz*tz)

	    ! Calculate integer bins
	    bin(1)=NINT(tx/binwidth)
	    bin(2)=NINT(ty/binwidth)
	    bin(3)=NINT(tz/binwidth)
	    bin(4)=NINT(dist/binwidth) - nbins

	    ! Store RMSD
	    sprmsdxyz(sp,1) = sprmsdxyz(sp,1) + tx*tx
	    sprmsdxyz(sp,2) = sprmsdxyz(sp,2) + ty*ty
	    sprmsdxyz(sp,3) = sprmsdxyz(sp,3) + tz*tz
	    sprmsdxyz(sp,4) = sprmsdxyz(sp,4) + dist*dist

	    ! Clamp bins to our bin range (so we can properly normalise the histograms)
	    do n=1,4
	      if (bin(n).gt.nbins) bin(n) = nbins
	      if (bin(n).lt.-nbins) bin(n) = -nbins
	      devxyz(sp,m,n,bin(n)) = devxyz(sp,m,n,bin(n)) + 1
	    end do

	    ! -----------------
	    ! Angular variation
	    ! -----------------

	    dp(1) = refaxes(sp,m,1)*axesA(sp,m,1) + refaxes(sp,m,2)*axesA(sp,m,2) + refaxes(sp,m,3)*axesA(sp,m,3)
	    dp(2) = refaxes(sp,m,4)*axesA(sp,m,4) + refaxes(sp,m,5)*axesA(sp,m,5) + refaxes(sp,m,6)*axesA(sp,m,6)
	    dp(3) = refaxes(sp,m,7)*axesA(sp,m,7) + refaxes(sp,m,8)*axesA(sp,m,8) + refaxes(sp,m,9)*axesA(sp,m,9)

	    ! Store RMSD
	    sprmsdangle(sp,1) = sprmsdangle(sp,1) + (1.0-dp(1))*(1.0-dp(1))
	    sprmsdangle(sp,2) = sprmsdangle(sp,2) + (1.0-dp(2))*(1.0-dp(2))
	    sprmsdangle(sp,3) = sprmsdangle(sp,3) + (1.0-dp(3))*(1.0-dp(3))
	    sprmsdangle(sp,4) = sprmsdangle(sp,4) + (1.0-dp(1))*(1.0-dp(1)) + (1.0-dp(2))*(1.0-dp(2)) + (1.0-dp(3))*(1.0-dp(3))

	    ! Clamp bins to bin range
	    do n=1,3
	      bin(n)=NINT(dp(n)/anglebinwidth)
	      if (bin(n).gt.nanglebins) bin(n) = nanglebins
	      if (bin(n).lt.-nanglebins) bin(n) = -nanglebins
	      devangle(sp,m,n,bin(n)) = devangle(sp,m,n,bin(n)) + 1
	    end do

	  end do
	end do

	! Next frame (or finish)
116	if (nframes.EQ.endf) goto 120
	goto 102

119	write(0,*) "Ended on file error."
120	write(0,*) "Finished."

	! Sum molecular quantities into species quantities
	do sp=1,nspecies
	  do m=1,s_nmols(sp)
	    spdevxyz(sp,1,:) = spdevxyz(sp,1,:) + devxyz(sp,m,1,:)
	    spdevxyz(sp,2,:) = spdevxyz(sp,2,:) + devxyz(sp,m,2,:)
	    spdevxyz(sp,3,:) = spdevxyz(sp,3,:) + devxyz(sp,m,3,:)
	    spdevxyz(sp,4,:) = spdevxyz(sp,4,:) + devxyz(sp,m,4,:)
	    spdevangle(sp,1,:) = spdevangle(sp,1,:) + devangle(sp,m,1,:)
	    spdevangle(sp,2,:) = spdevangle(sp,2,:) + devangle(sp,m,2,:)
	    spdevangle(sp,3,:) = spdevangle(sp,3,:) + devangle(sp,m,3,:)
	  end do
	end do

	! Normalise data
	do sp=1,nspecies
	  devxyz(sp,:,:,:) = devxyz(sp,:,:,:) / nframes
	  spdevxyz(sp,:,:) = spdevxyz(sp,:,:) / (nframes*s_nmols(sp))
	  sprmsdxyz(sp,:) = sqrt( sprmsdxyz(sp,:) / (nframes*s_nmols(sp)) )
	  devangle(sp,:,:,:) = devangle(sp,:,:,:) / nframes
	  spdevangle(sp,:,:) = spdevangle(sp,:,:) / (nframes*s_nmols(sp))
	  do n=1,3
	    sprmsdangle(sp,n) = sqrt( sprmsdangle(sp,n) / (nframes*s_nmols(sp)) )
	  end do
	  sprmsdangle(sp,4) = sqrt( sprmsdangle(sp,4) / ((nframes*s_nmols(sp))*3.0) )
	end do

	! Write out data
	do sp=1,nspecies
	  ! -----------------
	  ! Per molecule data
	  ! -----------------
	  do m=1,s_nmols(sp)
	    ! Construct destfile name
	    tth = m / 10000; n = m - tth*10000
	    th = n / 1000; n = n - th*1000
	    hun = n / 100; n = n - hun*100
	    ten = n / 10; n = n - ten*10
	    units = n
	    ! Positional variation
	    resfile=basename(1:baselen)//"devxyz_sp"//CHAR(48+sp)//"m"//CHAR(48+tth)//CHAR(48+th)//CHAR(48+hun)//CHAR(48+ten)//CHAR(48+units)
	    open(unit=11,file=resfile,form='formatted',status='replace')
	    write(11,*) "# DEVXYZ data for : species ",s_name(sp)
	    write(11,*) "#                   molecule ",m
	    write(11,*) "#   Bin   X    Y     Z"
	    do n=-nbins,nbins
	      write(11,"(4f10.3)") n*binwidth+0.5*binwidth, devxyz(sp,m,1,n), devxyz(sp,m,2,n), devxyz(sp,m,3,n)
	    end do
	    close(11)
	    ! Distance variation
	    resfile=basename(1:baselen)//"devdist_sp"//CHAR(48+sp)//"m"//CHAR(48+tth)//CHAR(48+th)//CHAR(48+hun)//CHAR(48+ten)//CHAR(48+units)
	    open(unit=11,file=resfile,form='formatted',status='replace')
	    write(11,*) "# DEVDIST data for : species ",s_name(sp)
	    write(11,*) "#                    molecule ",m
	    write(11,*) "#   Bin    Dist"
	    do n=-nbins,nbins
	      write(11,"(4f10.3)") (n+nbins)*binwidth+0.5*binwidth, devxyz(sp,m,4,n)
	    end do
	    close(11)
	    ! Angular variation
	    resfile=basename(1:baselen)//"devangle_sp"//CHAR(48+sp)//"m"//CHAR(48+tth)//CHAR(48+th)//CHAR(48+hun)//CHAR(48+ten)//CHAR(48+units)
	    open(unit=11,file=resfile,form='formatted',status='replace')
	    write(11,*) "# DEVANGLE data for : species ",s_name(sp)
	    write(11,*) "#                     molecule ",m
	    write(11,*) "#   Bin   X    Y     Z"
	    do n=-nanglebins,nanglebins
	      write(11,"(4f10.3)") n*anglebinwidth+0.5*anglebinwidth, devangle(sp,m,1,n), devangle(sp,m,2,n), devangle(sp,m,3,n)
	    end do
	    close(11)
	  end do
	  ! ----------------
	  ! Per species data
	  ! ----------------
	  ! Positional variation
	  resfile=basename(1:baselen)//"spdevxyz"//CHAR(48+sp)
	  open(unit=11,file=resfile,form='formatted',status='replace')
	  write(11,*) "# DEVXYZ data for species ",s_name(sp)
	  write(11,*) "#   Bin   X    Y     Z"
	  do n=-nbins,nbins
	    write(11,"(4f10.3)") n*binwidth+0.5*binwidth, spdevxyz(sp,1,n), spdevxyz(sp,2,n), spdevxyz(sp,3,n)
	  end do
	  close(11)
	  ! Distance variation
	  resfile=basename(1:baselen)//"spdevdist"//CHAR(48+sp)
	  open(unit=11,file=resfile,form='formatted',status='replace')
	  write(11,*) "# DEVDIST data for species ",s_name(sp)
	  write(11,*) "#   Bin    Dist"
	  do n=-nbins,nbins
	    write(11,"(4f10.3)") (n+nbins)*binwidth+0.5*binwidth, spdevxyz(sp,4,n)
	  end do
	  close(11)
	  ! Angle variation
	  resfile=basename(1:baselen)//"spdevangle"//CHAR(48+sp)
	  open(unit=11,file=resfile,form='formatted',status='replace')
	  write(11,*) "# DEVANGLE data for species ",s_name(sp)
	  write(11,*) "#   Bin    X    Y    Z"
	  do n=-nanglebins,nanglebins
	    write(11,"(4f10.3)") n*anglebinwidth+0.5*anglebinwidth, spdevangle(sp,1,n), spdevangle(sp,2,n), spdevangle(sp,3,n)
	  end do
	  close(11)
	  ! RMSD data
	  write(15,"(a,i2,a,4f10.5)") "Species ",sp,"     XYZR RMSD : ", (sprmsdxyz(sp,m),m=1,4)
	  write(15,"(a,i2,a,4f10.5)") "Species ",sp," XYZ(XYZ) RMSD : ", (sprmsdangle(sp,m),m=1,4)
	end do

	end program cryscomp
