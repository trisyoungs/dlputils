!	** glucanal **
!	Analyse the output from glucgeom

	program glucanal
	use dlprw; use utility
	implicit none
	character*80 :: comfile,hofile,flagfile
	character*20 :: temparg
	character*4 :: cutstr
	integer :: n,nframes,success,nargs,m,cl,gluc
	integer :: framestodo = 0, nsites, nchloride, nglucose, ndata, startf = -1
	integer :: hist_distbins, hist_anglebins, molcount, nmapframes, m1, m2
	integer :: iargc, norm
	logical :: molmap
	integer, allocatable :: cl_ownergluc(:,:), cl_glucsite(:,:)
	integer, allocatable :: molflags(:)
	real*8, allocatable :: cl_ohdist(:,:), cl_ohangle(:,:), cl_comdist(:,:)
	real*8, allocatable :: hist_distcom(:), hist_ohdist1(:), hist_ohangle(:), hist_ohdist2(:), hist_ohdist5(:)
	real*8, allocatable :: hist_2ndohdist(:), hist_2ndohangle(:), hist_2ndohdist_cut(:), hist_2ndohangle_cut(:)
	real*8, allocatable :: hist_ohgroup(:,:), hist_ohgrouppairs(:,:)
	real*8 :: hist_distbin, hist_anglebin, hist_distmax, hist_anglemax, ohcut=20.0

	nargs = iargc()
	if (nargs.LT.4) stop "Usage : glucanal <com.dat> <ho.dat> <nchloride> <ngluc> [-cutoff cut] [-framestodo n] [-start n]"
	call getarg(1,comfile)
	call getarg(2,hofile)
	call getarg(3,temparg)
	read(temparg,"(I5)") nchloride
	call getarg(4,temparg)
	read(temparg,"(I5)") nglucose
	
	n = 4
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temparg)
	  select case (temparg)
	    case ("-cutoff")
	      n = n + 1; call getarg(n,temparg); read(temparg,"(f20.10)") ohcut
	    case ("-frames")
	      n = n + 1; call getarg(n,temparg); read(temparg,"(I6)") framestodo
	      write(0,"(A,I6)") "Number of frames to sample: ",framestodo
	    case ("-start")
	      n = n + 1; call getarg(n,temparg); read(temparg,"(I6)") startf
	      write(0,"(A,I4)") "Starting frame = ",startf
	    case ("-molmap")
	      n = n + 1; call getarg(n,flagfile)
	      write(0,"(A,A)") "Using flags for species 1 molecules from : ",flagfile
	      open(unit=25,file=flagfile,form="formatted",status="old")
	      molmap = .TRUE.
	    case default
	      write(0,*) "Unrecognised command line option : ",temparg
	      stop
	  end select
	end do

	write(0,"(A,F5.2,A)") "Chloride-OH interaction cutoff is ",ohcut," Angstroms"

	! Open
	open(unit=19,file=comfile,form="formatted",status="old")
	open(unit=20,file=hofile,form="formatted",status="old")

	! Set some numbers
	nsites = 5
	if (nglucose.eq.1) then
	  ndata = 5
	else
	  ndata = 10
	end if
	hist_distbin = 0.05
	hist_distmax = 15.0
	hist_anglebin = 0.5
	hist_anglemax = 180.0
	hist_distbins = hist_distmax / hist_distbin
	hist_anglebins = hist_anglemax / hist_anglebin

	if (molmap) allocate(molflags(nglucose))
	allocate(cl_comdist(nchloride,nglucose))
	allocate(cl_ohdist(nchloride,ndata))
	allocate(cl_ownergluc(nchloride,ndata))
	allocate(cl_ohangle(nchloride,ndata))
	allocate(cl_glucsite(nchloride,ndata))

	allocate(hist_distcom(hist_distbins))
	allocate(hist_ohdist1(hist_distbins))
	allocate(hist_ohdist2(hist_distbins))
	allocate(hist_ohdist5(hist_distbins))
	allocate(hist_ohangle(hist_anglebins))
	allocate(hist_2ndohdist(hist_distbins))
	allocate(hist_2ndohangle(hist_anglebins))
	allocate(hist_2ndohdist_cut(hist_distbins))
	allocate(hist_2ndohangle_cut(hist_anglebins))

	allocate(hist_ohgroup(3,5))
	allocate(hist_ohgrouppairs(5,5))

	! Initialise variables
	hist_distcom = 0.0
	hist_ohdist1 = 0.0
	hist_ohdist2 = 0.0
	hist_ohdist5 = 0.0
	hist_ohangle = 0.0
	hist_2ndohdist = 0.0
	hist_2ndohangle = 0.0
	hist_ohgroup = 0.0
	hist_ohgrouppairs = 0.0

	! XXXX
	! XXXX Main routine....
	! XXXX
	nframes=0
	nmapframes=0
	! If we're using a mapfile, read it here
101     if (molmap) then
	  read(25,"(I5,1x,40(I2,1x))",end=118,err=117) m1,(molflags(m2),m2=1,nglucose)
	  ! if (m1.NE.nframes) stop "Mapfile missed a frame number!"
	end if
102	read(19,"(6x,A20)",end=801,err=799) temparg
	read(20,"(6x,A20)",end=801,err=799) temparg
	do n=1,nchloride
	  ! Read in com data
	  read(19,"(4x,20f7.3)",end=801,err=799) (cl_comdist(n,m),m=1,nglucose)
	  ! Read in oh data
	  read(20,"(4x,10(i2,1x,i1,1x,f7.3,f7.2,1x))",end=801,err=799) (cl_ownergluc(n,m),cl_glucsite(n,m),cl_ohdist(n,m),cl_ohangle(n,m),m=1,ndata)
	end do
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) then
	  if (molmap) then
	    write(0,"(I6,2x,'(',I6,',',I6,')')") nframes,nmapframes,molcount
	  else
	    write(0,"(i6)") nframes
	  end if
	end if

	if (nframes.LT.startf) goto 102

	if ((molmap).and.(m1.ne.nframes)) goto 102
	if (molmap) nmapframes = nmapframes + 1

	! ANALYSIS STARTS HERE

	! Histogram of Cl-gluc[com] distances
	do gluc=1,nglucose
	  if (molmap.and.(molflags(gluc).eq.0)) cycle
	  do cl=1,nchloride
	    molcount = molcount + 1
	    m = (cl_comdist(cl,gluc) / hist_distbin) + 1
	    if (m.lt.hist_distbins) hist_distcom(m) = hist_distcom(m) + 1
	  end do
	end do

	do cl=1,nchloride
	  ! Histograms of Cl-OH distances below ohcut, and the corresponding angle distribution
	  ! Consider best contact, best two contacts etc.
	  do n=1,ndata
	    if (molmap.and.(molflags(cl_ownergluc(cl,n)).eq.0)) cycle
	    if (cl_ohdist(cl,n).le.ohcut) then
	      ! Bin Cl-OH distance
	      m = (cl_ohdist(cl,n) / hist_distbin) + 1
	      if (m.lt.hist_distbins) then
		if (n.eq.1) then
		  hist_ohdist1(m) = hist_ohdist1(m) + 1
		  hist_ohdist2(m) = hist_ohdist2(m) + 1
		  hist_ohdist5(m) = hist_ohdist5(m) + 1
		end if
		if (n.gt.1) then
		  hist_ohdist2(m) = hist_ohdist2(m) + 1
		  hist_ohdist5(m) = hist_ohdist5(m) + 1
		end if
		if (n.gt.2) hist_ohdist5(m) = hist_ohdist5(m) + 1
	      end if
	      ! Bin Cl-OH angle
	      m = (cl_ohangle(cl,n) / hist_anglebin) + 1
	      if (m.lt.hist_anglebins) hist_ohangle(m) = hist_ohangle(m) + 1
	    end if
	  end do
	  ! Histogram of second bond distances and angles where the shortest contact is less than ohcut
	  if (cl_ohdist(cl,1).le.ohcut) then
	    if (molmap.and.(molflags(cl_ownergluc(cl,1)).ne.0)) then
	      m = (cl_ohdist(cl,2) / hist_distbin) + 1
	      n = (cl_ohangle(cl,2) / hist_anglebin) + 1
	      ! Always add to general histograms
	      if (m.lt.hist_distbins) then
	        hist_2ndohdist(m) = hist_2ndohdist(m) + 1
	        hist_2ndohangle(n) = hist_2ndohangle(n) + 1
	      end if
	      if (cl_ohdist(cl,2).le.ohcut) then
	        hist_2ndohdist_cut(m) = hist_2ndohdist_cut(m) + 1
	        hist_2ndohangle_cut(n) = hist_2ndohangle_cut(n) + 1
	      end if
	    end if
	  end if
	  ! Histogram of OH groups involved in bonds with chloride, below the cutoff distance
	  do n=1,ndata
	    if (cl_ohdist(cl,n).le.ohcut) then
	      if (molmap.and.(molflags(cl_ownergluc(cl,n)).eq.0)) cycle
	      m = cl_glucsite(cl,n)
	      if (n.eq.1) then
		hist_ohgroup(1,m) = hist_ohgroup(1,m) + 1
		hist_ohgroup(2,m) = hist_ohgroup(2,m) + 1
		hist_ohgroup(3,m) = hist_ohgroup(3,m) + 1
	      end if
	      if (n.gt.1) then
		hist_ohgroup(2,m) = hist_ohgroup(2,m) + 1
		hist_ohgroup(3,m) = hist_ohgroup(3,m) + 1
	      end if
	      if (n.gt.2) hist_ohgroup(3,m) = hist_ohgroup(3,m) + 1
	    end if
	  end do
	  ! Correlation between OH group pairs for TWO chlorides less than ohcut angleance away
	  if (cl_ohdist(cl,1).le.ohcut) then
	    if (cl_ohdist(cl,2).le.ohcut) then
	      if (molmap) then
		if (molflags(cl_ownergluc(cl,1)).eq.0) cycle
		if (molflags(cl_ownergluc(cl,2)).eq.0) cycle
	      end if
	      m = cl_glucsite(cl,1)
	      n = cl_glucsite(cl,2)
		! HERE, distinguish between OH groups on same / different glucose molecules
	      hist_ohgrouppairs(m,n) = hist_ohgrouppairs(m,n) + 1
	    end if
	  end if
	end do


	if (nframes.EQ.framestodo) goto 800
	! Next frame
	goto 101

117     write(0,*) "Error reading mapfile."
        write(0,*) "Selected ",nmapframes," from trajectory before error."
        goto 801
118     write(0,*) "Reached end of mapfile."
        write(0,*) "Selected ",nmapframes," from trajectory."
        goto 801
799	write(0,*) "Error in data files."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",framestodo,")"
	goto 801
800	write(0,*) "Framestodo was fulfilled."
801	write(0,*) ""

	! WRITE OUTPUT

	! Histogram of Cl-gluc[com] distances
	open(unit=14,file="glucanal_cl-gluccom.dat",form="formatted",status="replace")
	do n=1,hist_distbins
	  write(14,"(f7.3,2x,f10.6)") hist_distbin*n-0.5*hist_distbin,hist_distcom(n) / nframes
	end do
	close(14)

	if (molmap) then
	  norm = nmapframes
	else
	  norm = nframes
	end if
	! Histogram of Cl-OH distances below ohcut, and the corresponding angle distribution
	open(unit=14,file="glucanal_cl-oh1-"//cutstr//".dat",form="formatted",status="replace")
	open(unit=15,file="glucanal_cl-oh2-"//cutstr//".dat",form="formatted",status="replace")
	open(unit=16,file="glucanal_cl-oh5-"//cutstr//".dat",form="formatted",status="replace")
	open(unit=17,file="glucanal_cl-2ndoh.dat",form="formatted",status="replace")
	open(unit=18,file="glucanal_cl-2ndoh-"//cutstr//".dat",form="formatted",status="replace")
	do n=1,hist_distbins
	  write(14,"(f7.3,2x,f10.6)") hist_distbin*n-0.5*hist_distbin,hist_ohdist1(n) / norm
	  write(15,"(f7.3,2x,f10.6)") hist_distbin*n-0.5*hist_distbin,hist_ohdist2(n) / norm
	  write(16,"(f7.3,2x,f10.6)") hist_distbin*n-0.5*hist_distbin,hist_ohdist5(n) / norm
	  write(17,"(f7.3,2x,f10.6)") hist_distbin*n-0.5*hist_distbin,hist_2ndohdist(n) / norm
	  write(18,"(f7.3,2x,f10.6)") hist_distbin*n-0.5*hist_distbin,hist_2ndohdist_cut(n) / norm
	end do
	close(14)
	close(15)
	close(16)
	close(17)
	close(18)

	open(unit=14,file="glucanal_cl-ohangle-"//cutstr//".dat",form="formatted",status="replace")
	open(unit=15,file="glucanal_cl-2ndohangle.dat",form="formatted",status="replace")
	open(unit=16,file="glucanal_cl-2ndohangle-"//cutstr//".dat",form="formatted",status="replace")
	do n=1,hist_anglebins
	  write(14,"(f7.3,2x,f10.6)") hist_anglebin*n-0.5*hist_anglebin,hist_ohangle(n) / norm
	  write(15,"(f7.3,2x,f10.6)") hist_anglebin*n-0.5*hist_anglebin,hist_2ndohangle(n) / norm
	  write(16,"(f7.3,2x,f10.6)") hist_anglebin*n-0.5*hist_anglebin,hist_2ndohangle_cut(n) / norm
	end do
	close(14)
	close(15)
	close(16)

	open(unit=14,file="glucanal_cl-ohgroup-"//cutstr//".dat",form="formatted",status="replace")
	open(unit=15,file="glucanal_cl-ohgroupcorr-"//cutstr//".dat",form="formatted",status="replace")
	write(15,"(9x,5i10)") (n,n=1,5)
	do n=1,5
	  write(14,"(i2,2x,5f10.6)") n,(hist_ohgroup(m,n) / nframes,m=1,3)
	  write(15,"(i2,2x,5f10.6)") n,(hist_ohgrouppairs(m,n) / nframes,m=1,5)
	end do
	close(14)

999	write(0,*) "Finished."
	end program glucanal


