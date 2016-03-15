!	** bident4 **
!	Calculate style of binding between 'bidentate' species and positions on outer species

	program bident4
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,headerfile, strtmp
	character*20 :: temparg
	integer, parameter :: MAXSITES = 20
	integer :: MAXATOMCONTACTS = 10, MAXPATTERNS = 30000, MAXMOLCONTACTS = 6
	integer :: nframes,success,nargs,n,m1,m2,i,j,k,baselen,aoff1,aoff2,m
	integer :: framestodo = -1, framestodiscard = 0, sp1, sp2, nsp2sites(MAXSP), sp2sites(MAXSP,MAXSITES,2), sp1sites(2)
	integer :: iargc, site, site2, closebit, bit, npatterns = 0, ncontacts, found, t1, t2, t3
	integer :: tth,th,hun,ten,units,natoms, grid = 40, binx, biny, binz, mtemp(3)
	integer, allocatable :: centralbits(:), nsp1molcontacts(:), sp1molcontacts(:,:,:), nsp1atomcontacts(:), sp2molcontacts(:,:), nsp2molcontacts(:)
	integer :: nsp2atomcontacts(MAXSP,MAXSITES,2), writecontact = -1
	logical :: done, calc3d = .FALSE., writeone = .FALSE., siteequiv = .FALSE.
	real*8 :: a(3),b(3),c(3),mima(3),mimb(3),dac,dbc,total,avgmol,avgatom,tx,ty,tz,px,py,pz
	real*8 :: delta = 0.25, sp2sitemaxdist(MAXSP,MAXSITES), cp
	integer, allocatable :: pattern(:,:,:), nmolcontacts(:), natomcontacts(:), pfreq(:), ptemp(:,:)
	real*8, allocatable :: geom(:,:,:), gtemp(:,:), pdens(:,:,:,:), pdenstemp(:,:,:)

	nargs = iargc()
	if (nargs.lt.6) then
	  write(0,"(a120)") "Usage : bident4 <HISTORYfile> <OUTfile> <sp1> <sp1 atom i> <sp1 atom j>"
	  write(0,"(a80)") "    [-frames n] [-discard n] [-site sp2 a b maxdist (e.g. H(a)-O(b))] [-maxpatterns n]"
	  write(0,"(a80)") "    [-grid n] [-delta d] [-maxatomcontacts n] [-maxmolcontacts n]"
	  write(0,"(a80)") "    [-calc3d x1 x2 y1 y2] [-writeone] [-writecontact sp]"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temparg); read(temparg,"(I5)") sp1
	call getarg(4,temparg); read(temparg,"(I5)") sp1sites(1)
	call getarg(5,temparg); read(temparg,"(I5)") sp1sites(2)
	
	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	call alloc_axis()
	nsp2sites = 0
	sp2sites = 0

	n = 5
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temparg)
          select case (temparg)
            case ("-site")
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") sp2
	      nsp2sites(sp2) = nsp2sites(sp2) + 1
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") sp2sites(sp2,nsp2sites(sp2),1)
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") sp2sites(sp2,nsp2sites(sp2),2)
              n = n + 1; call getarg(n,temparg); read(temparg,"(f12.6)") sp2sitemaxdist(sp2,nsp2sites(sp2))
              write(0,"(A,i2,A,i4,' and ',i4,a,f12.6)") "Site on species ",sp2," defined: atoms are ",(sp2sites(sp2,nsp2sites(sp2),m),m=1,2), "dist = ", sp2sitemaxdist(sp2,nsp2sites(sp2))
            case ("-header")
              n = n + 1; call getarg(n,headerfile)
            case ("-frames")
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") framestodo
            case ("-discard")
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") framestodiscard
            case ("-maxatomcontacts")
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") MAXATOMCONTACTS
	      write(0,*) "MAXATOMCONTACTS set to ", MAXATOMCONTACTS
            case ("-maxmolcontacts")
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") MAXMOLCONTACTS
	      write(0,*) "MAXMOLCONTACTS set to ", MAXMOLCONTACTS
            case ("-maxpatterns")
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") MAXPATTERNS
	      write(0,*) "MAXPATTERNS set to ", MAXPATTERNS
	    case ("-calc3d") 
	      n = n + 1; call getarg(n,temparg); read(temparg,"(I3)") axesAatoms(sp1,1)
	      n = n + 1; call getarg(n,temparg); read(temparg,"(I3)") axesAatoms(sp1,2)
	      n = n + 1; call getarg(n,temparg); read(temparg,"(I3)") axesAatoms(sp1,3)
	      n = n + 1; call getarg(n,temparg); read(temparg,"(I3)") axesAatoms(sp1,4)
	      write(0,"(A,I2,A,I2,A,I2,A,I2,A)") "Local axes for central species calculated from: X=",axesAatoms(sp1,1),"->", &
	        & axesAatoms(sp1,2),", Y=0.5X->0.5(r(",axesAatoms(sp1,3),")->r(",axesAatoms(sp1,4),"))"
	      axesAdefined(sp1) = .true.
	      calc3d = .TRUE.
	    case ("-equiv")
	      write(0,*) "Sites on central species will be treated as equivalent for single contacts."
	      siteequiv = .TRUE.
	    case ("-writeone")
	      writeone = .TRUE.
	    case ("-writecontact")
	      n = n + 1; call getarg(n,temparg); read(temparg,"(I3)") writecontact 
	      write(0,*) "Contact patterns around central will be printed per molecule, per frame", writecontact
            case ("-grid")
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") grid
	      write(0,*) "Grid points in each +/- direction : ", grid
            case ("-delta")
              n = n + 1; call getarg(n,temparg); read(temparg,"(f12.5)") delta
	      write(0,*) "Grid point spacing : ", delta
	    case default
	      write(0,*) "Unrecognised argument:", temparg
	      stop
	  end select
	end do
	
	! Open output files
	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="rdfresults."
	  baselen=11
	else
	  basename=hisfile(1:baselen)
	endif
	open(unit=11,file=basename(1:baselen)//"bident4_"//CHAR(48+sp1),form='formatted',status='replace')
	write(11,"('# Central species is ',i4)") sp1
	write(11,"('# Atom IDs used for interaction points :',i4,i4)") sp1sites
	do sp2=1,nspecies
	  write(11,"('# Number of sites defined on species ',i4)") nsp2sites(sp2)
	  do n=1,nsp2sites(sp2)
	    write(11,"(6x,2i5)") sp2sites(sp2,n,1),sp2sites(sp2,n,2)
	  end do
	end do

	write(11,*) "Max number of contacts is ", maxval(nsp2sites)

	allocate(centralbits(2**(maxval(nsp2sites)*2)))
	allocate(sp1molcontacts(s_nmols(sp1),MAXMOLCONTACTS,3))
	allocate(nsp1atomcontacts(s_nmols(sp1)))
	allocate(nsp1molcontacts(s_nmols(sp1)))
	if (writecontact.gt.0) then
	  allocate(sp2molcontacts(s_nmols(writecontact),MAXMOLCONTACTS))
	  allocate(nsp2molcontacts(s_nmols(writecontact)))
	end if
	allocate(pattern(MAXPATTERNS,MAXMOLCONTACTS,3))
	allocate(ptemp(MAXMOLCONTACTS,3))
	allocate(pfreq(MAXPATTERNS))
	allocate(nmolcontacts(MAXPATTERNS))
	allocate(natomcontacts(MAXPATTERNS))
	if (calc3d) then
	  allocate(geom(MAXPATTERNS, s_natoms(sp1)+maxval(s_natoms)*MAXMOLCONTACTS, 3), gtemp(s_natoms(sp1)+maxval(s_natoms)*MAXMOLCONTACTS, 3))
	  allocate(pdens(MAXPATTERNS, -grid:grid, -grid:grid, -grid:grid), pdenstemp(-grid:grid, -grid:grid, -grid:grid)) 
	end if

	pattern = 0
	pfreq = 0
	nmolcontacts = 0
	natomcontacts = 0
	nsp2atomcontacts = 0	
	if (calc3d) then
	  geom = 0.0
	  pdens = 0.0
	end if

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	if (framestodiscard.ne.0) then
	  framestodiscard = framestodiscard - 1
	  goto 101
	end if
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	! Generate all molecular axes
	if (calc3d) call genaxes()

	! Loop over central species (bidentate target)

	nsp1atomcontacts = 0
	nsp1molcontacts = 0
	sp1molcontacts = 0
	if (writecontact.gt.0) then
	  nsp2molcontacts = 0
	  sp2molcontacts = 0
	end if

	aoff1 = s_start(sp1)-1
	do m1 = 1,s_nmols(sp1)

	  ! Grab coordinates of central species site pair
	  a(1)=xpos(aoff1 + sp1sites(1))
	  a(2)=ypos(aoff1 + sp1sites(1))
	  a(3)=zpos(aoff1 + sp1sites(1))
	  b(1)=xpos(aoff1 + sp1sites(2))
	  b(2)=ypos(aoff1 + sp1sites(2))
	  b(3)=zpos(aoff1 + sp1sites(2))

	  ! Now, loop over all species, skipping the central one, searching for contacts
	  do sp2=1,nspecies

	    if (sp1.eq.sp2) cycle

	    aoff2 = s_start(sp2)-1
	    do m2=1,s_nmols(sp2)

	      closebit = 0
	      ncontacts = 0

	      ! For each outer molecule, go over all defined sites on it and encode a bitmask of the contacts
	      do site = 1,nsp2sites(sp2)

	        ! Grab this central site position
	        c(1)=xpos(aoff2+sp2sites(sp2,site,1))
	        c(2)=ypos(aoff2+sp2sites(sp2,site,1))
	        c(3)=zpos(aoff2+sp2sites(sp2,site,1))

	        ! Get MIM positions of outer species sites w.r.t. this site
	        call pbc(a(1),a(2),a(3),c(1),c(2),c(3),mima(1),mima(2),mima(3))
	        call pbc(b(1),b(2),b(3),c(1),c(2),c(3),mimb(1),mimb(2),mimb(3))
	        mima = mima - c
	        mimb = mimb - c

	        ! Determine distances
	        dac = sqrt(mima(1)*mima(1) + mima(2)*mima(2) + mima(3)*mima(3))
	        dbc = sqrt(mimb(1)*mimb(1) + mimb(2)*mimb(2) + mimb(3)*mimb(3))

	        ! Calculate closeness bit for central molecule
		bit = 0
	        if (dac.lt.sp2sitemaxdist(sp2,site)) then
		  bit = bit + 2**(site-1)
		  nsp2atomcontacts(sp2,site,1) = nsp2atomcontacts(sp2,site,1) + 1
		  ncontacts = ncontacts + 1
		end if
		if ((sp1sites(1).ne.sp1sites(2)).and.(dbc.lt.sp2sitemaxdist(sp2,site))) then
		  bit = bit + 2**(nsp2sites(sp2)+site-1)
		  nsp2atomcontacts(sp2,site,2) = nsp2atomcontacts(sp2,site,2) + 1
		  ncontacts = ncontacts + 1
		end if

		closebit = closebit + bit

	      end do

 
	      ! Equivalency of sites in single and bifurcated contacts
	      if (siteequiv) then
	        if (ncontacts.eq.1) then
	          ! Single contact - use lower of the two possible bits
	          if (closebit.gt.(2**(nsp2sites(sp2)-1))) closebit = closebit / (2**nsp2sites(sp2))
	        else if (ncontacts.eq.2) then
	          ! Bidentate interactions need no adjustment
		  done = .FALSE.
		  do site = 0,nsp2sites(sp2)-2
		    do site2 = site+1,nsp2sites(sp2)-1
		      ! Bridging interactions are equivalent whichever way round they are... choose smaller number...
		      if (closebit.eq.(2**(site)+2**(site2+nsp2sites(sp2)))) then
			!write(0,*) "OLD,NEW",closebit,2**site2 + 2**(site+nsp2sites(sp2))
			closebit = 2**site2 + 2**(site+nsp2sites(sp2))
			done = .TRUE.
			exit
		      end if
		      ! Bifurcated interactions can use the smaller of the two possible values as well
		      if (closebit.eq.(2**(site+nsp2sites(sp2))+2**(site2+nsp2sites(sp2)))) then
			!write(0,*) "BIFOLD,NEW",closebit,2**site2 + 2**site
			closebit = 2**site + 2**site2
			done = .TRUE.
			exit
		      end if
		    end do
		    if (done) exit
		  end do
	          ! ???
	        end if
	      end if

	      ! Does this central molecule have any contacts with this outer molecule?
	      if (ncontacts.gt.0) then
		nsp1molcontacts(m1) = nsp1molcontacts(m1)+1
		if (nsp1molcontacts(m1).gt.MAXMOLCONTACTS) stop "MAXMOLCONTACTS exceeded."
		nsp1atomcontacts(m1) = nsp1atomcontacts(m1)+ncontacts
	        sp1molcontacts(m1,nsp1molcontacts(m1),1) = sp2
	        sp1molcontacts(m1,nsp1molcontacts(m1),2) = m2
	        sp1molcontacts(m1,nsp1molcontacts(m1),3) = closebit
	        ! Store contact info around central molecule
		if (writecontact.eq.sp2) then
		  nsp2molcontacts(m2) = nsp2molcontacts(m2)+1
		  if (nsp2molcontacts(m2).gt.MAXMOLCONTACTS) stop "MAXMOLCONTACTS exceeded."
	          sp2molcontacts(m2,nsp2molcontacts(m2)) = closebit
		end if
	      end if

	      aoff2 = aoff2 + s_natoms(sp2)

	    end do

	  end do

!	write(0,*) "Found ", nsp1molcontacts(m1), "for central molecule ", m1

	  aoff1 = aoff1 + s_natoms(sp1)

	end do

	! Write out sp2 contact pattern
	if (writecontact.gt.0) then
	  do m2=1,s_nmols(writecontact)
	    write(55,"(20i6)") m2,sp2molcontacts(m2,1:nsp2molcontacts(m2))
	  end do
	end if

	! Enumerate data
	do m1=1,s_nmols(sp1)

	  ! Sort contact data by key data
!	  write(0,"(a,50i4)") "Data was,",(sp1molcontacts(m1,n,1),sp1molcontacts(m1,n,2),sp1molcontacts(m1,n,3),n=1,nsp1molcontacts(m1))
	  do sp2 = 1,nspecies
	    ! Find first and last molecular contact involving this species
	    do i=1,nsp1molcontacts(m1)
	      if (sp1molcontacts(m1,i,1).eq.sp2) exit
	    end do
	    if (i.gt.nsp1molcontacts(m1)) cycle
	    do j=i,nsp1molcontacts(m1)
	      if (sp1molcontacts(m1,j,1).ne.sp2) exit
	    end do
	    j = j - 1
	    ! Bubble sort
	    do n=i,j
	      do m=i,j-1
	        if (sp1molcontacts(m1,m+1,3).lt.sp1molcontacts(m1,m,3)) then
		  mtemp(:) = sp1molcontacts(m1,m,:)
		  sp1molcontacts(m1,m,:) = sp1molcontacts(m1,m+1,:)
		  sp1molcontacts(m1,m+1,:) = mtemp(:)
		end if
	      end do
	    end do
	  end do
!	  write(0,"(a,50i4)") "------> ",(sp1molcontacts(m1,n,1),sp1molcontacts(m1,n,2),sp1molcontacts(m1,n,3),n=1,nsp1molcontacts(m1))

	  ! Each central molecule now has a contact pattern. Search for each to see if we've encountered it before.
	  ! If not, add a new pattern into the list

	  found = 0
! 	  write(0,"(a,50i4)") "Contact data is,",(sp1molcontacts(m1,n,1),sp1molcontacts(m1,n,2),sp1molcontacts(m1,n,3),n=1,nsp1molcontacts(m1))
	  do n=1,npatterns
	    if (nmolcontacts(n).ne.nsp1molcontacts(m1)) cycle
	    if (natomcontacts(n).ne.nsp1atomcontacts(m1)) cycle
	    !if (nmolcontacts(n).eq.0) found = n
	    t1 = 0
	    do i=1,nmolcontacts(n)
	      if (sp1molcontacts(m1,i,1).ne.pattern(n,i,1)) cycle
	      !if (sp1molcontacts(m1,i,2).ne.pattern(n,i,2)) cycle
	      if (sp1molcontacts(m1,i,3).ne.pattern(n,i,3)) cycle
	      t1 = t1 + 1
	    end do
	    ! Did we succesfully match all contacts?
	    if (t1.eq.nmolcontacts(n)) then
	      found = n
	      exit
	    end if
	  end do

	  if (found.ne.0) then
!  	    write(0,"(a,50i4)") "==== Matched   ,",(pattern(found,i,1),pattern(found,i,2),pattern(found,i,3),i=1,nmolcontacts(found))
	    pfreq(found) = pfreq(found) + 1
	  else
! 	    write(0,"(a)") "!No match"
	    ! Create new pattern
	    npatterns = npatterns + 1
	    if (npatterns.gt.MAXPATTERNS) stop "MAXPATTERNS exceeded. Adjust and run again!"
	    pfreq(npatterns) = 1
	    nmolcontacts(npatterns) = nsp1molcontacts(m1)
	    natomcontacts(npatterns) = nsp1atomcontacts(m1)
	    do i=1,nsp1molcontacts(m1)
	      pattern(npatterns,i,1) = sp1molcontacts(m1,i,1)
	      pattern(npatterns,i,2) = sp1molcontacts(m1,i,2)
	      pattern(npatterns,i,3) = sp1molcontacts(m1,i,3)
	    end do
	    found = npatterns
	  end if

	  ! Update 3D geometry if requested
	  if (calc3d.and.((writeone.and.(pfreq(found).eq.1)).or.(.not.writeone))) then
	    ! Central molecule is first
! 	write(0,*) "adding geometry to pattern geom() ",found
	    aoff1 = s_start(sp1) + s_natoms(sp1)*(m1-1) - 1
	!write(0,"(4f12.6)") axisx(sp1,m1,1), axisx(sp1,m1,2), axisx(sp1,m1,3), sqrt(sum(axisx(sp1,m1,:)*axisx(sp1,m1,:)))
	!write(0,"(4f12.6)") axisy(sp1,m1,1), axisy(sp1,m1,2), axisy(sp1,m1,3), sqrt(sum(axisy(sp1,m1,:)*axisy(sp1,m1,:)))
	!write(0,"(4f12.6)") axisz(sp1,m1,1), axisz(sp1,m1,2), axisz(sp1,m1,3), sqrt(sum(axisz(sp1,m1,:)*axisz(sp1,m1,:)))
	    do i=1,s_natoms(sp1)
	      ! Calculate minimum image position of each atom about axis origin
	      call pbc(xpos(aoff1+i),ypos(aoff1+i),zpos(aoff1+i),axesAorigin(sp1,m1,1),axesAorigin(sp1,m1,2),axesAorigin(sp1,m1,3),tx,ty,tz)
	      px = tx - axesAorigin(sp1,m1,1)
	      py = ty - axesAorigin(sp1,m1,2)
	      pz = tz - axesAorigin(sp1,m1,3)
	      tx = px*axesA(sp1,m1,1) + py*axesA(sp1,m1,2) + pz*axesA(sp1,m1,3)
	      ty = px*axesA(sp1,m1,4) + py*axesA(sp1,m1,5) + pz*axesA(sp1,m1,6)
	      tz = px*axesA(sp1,m1,7) + py*axesA(sp1,m1,8) + pz*axesA(sp1,m1,9)
	      geom(found,i,1) = geom(found,i,1) + tx
	      geom(found,i,2) = geom(found,i,2) + ty
	      geom(found,i,3) = geom(found,i,3) + tz
	    end do
	    ! Now surrounding molecules
	    aoff1 = s_natoms(sp1)
	    do n=1,nsp1molcontacts(m1)
	      sp2 = sp1molcontacts(m1,n,1)
	      m2 = sp1molcontacts(m1,n,2)
	      aoff2 = s_start(sp2) + s_natoms(sp2)*(m2-1) - 1
	      do i = 1,s_natoms(sp2)
		! Calculate minimum image position of each atom about central molecule axis origin
		call pbc(xpos(aoff2+i),ypos(aoff2+i),zpos(aoff2+i),axesAorigin(sp1,m1,1),axesAorigin(sp1,m1,2),axesAorigin(sp1,m1,3),tx,ty,tz)
		px = tx - axesAorigin(sp1,m1,1)
		py = ty - axesAorigin(sp1,m1,2)
		pz = tz - axesAorigin(sp1,m1,3)
		tx = px*axesA(sp1,m1,1) + py*axesA(sp1,m1,2) + pz*axesA(sp1,m1,3)
		ty = px*axesA(sp1,m1,4) + py*axesA(sp1,m1,5) + pz*axesA(sp1,m1,6)
		tz = px*axesA(sp1,m1,7) + py*axesA(sp1,m1,8) + pz*axesA(sp1,m1,9)
		geom(found,i+aoff1,1) = geom(found,i+aoff1,1) + tx
		geom(found,i+aoff1,2) = geom(found,i+aoff1,2) + ty
		geom(found,i+aoff1,3) = geom(found,i+aoff1,3) + tz
	      end do
	      ! Update pdens grid
	      do site=1,nsp2sites(sp2)
		if (btest(sp1molcontacts(m1,n,3),site-1).or.btest(sp1molcontacts(m1,n,3),nsp2sites(sp2)+site-1)) then
		  i = sp2sites(sp2,site,1)
		  call pbc(xpos(aoff2+i),ypos(aoff2+i),zpos(aoff2+i),axesAorigin(sp1,m1,1),axesAorigin(sp1,m1,2),axesAorigin(sp1,m1,3),tx,ty,tz)
		  px = tx - axesAorigin(sp1,m1,1)
		  py = ty - axesAorigin(sp1,m1,2)
		  pz = tz - axesAorigin(sp1,m1,3)
		  tx = px*axesA(sp1,m1,1) + py*axesA(sp1,m1,2) + pz*axesA(sp1,m1,3)
		  ty = px*axesA(sp1,m1,4) + py*axesA(sp1,m1,5) + pz*axesA(sp1,m1,6)
		  tz = px*axesA(sp1,m1,7) + py*axesA(sp1,m1,8) + pz*axesA(sp1,m1,9)
		  binx = NINT(tx/delta)
		  biny = NINT(ty/delta)
		  binz = NINT(tz/delta)
		  if ((abs(binx).gt.grid).or.(abs(biny).gt.grid).or.(abs(binz).gt.grid)) stop "Increase grid size and rerun."
		  pdens(found,binx,biny,binz) = pdens(found,binx,biny,binz) + 1.0
		end if
	      end do
	      aoff1 = aoff1 + s_natoms(sp2)
	    end do
	  end if

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

	! Sort pattern data by frequency (via Bubble sort)
	do n=1,npatterns
	  do m=1,npatterns-1
	    if (pfreq(m+1).gt.pfreq(m)) then
	      ! Swap data in 'm' with data in 'm+1'
	      ptemp(:,:) = pattern(m,:,:)
	      t1 = pfreq(m)
	      t2 = nmolcontacts(m)
	      t3 = natomcontacts(m)
	      pattern(m,:,:) = pattern(m+1,:,:)
	      pfreq(m) = pfreq(m+1)
	      nmolcontacts(m) = nmolcontacts(m+1)
	      natomcontacts(m) = natomcontacts(m+1)
	      pattern(m+1,:,:) = ptemp(:,:)
	      pfreq(m+1) = t1
	      nmolcontacts(m+1) = t2
	      natomcontacts(m+1) = t3
	      gtemp(:,:) = geom(m,:,:)
	      geom(m,:,:) = geom(m+1,:,:)
	      geom(m+1,:,:) = gtemp(:,:)
	      pdenstemp(:,:,:) = pdens(m,:,:,:)
	      pdens(m,:,:,:) = pdens(m+1,:,:,:)
	      pdens(m+1,:,:,:) = pdenstemp(:,:,:)
	    end if
	  end do
	end do

	! Output key data
	write(11,*) ""
	write(11,*) "*** Key data for bitmasks"
	write(11,*) ""
	do sp2=1,nspecies
	  if (sp1.eq.sp2) cycle
	  write(11,*) "Species ",sp2
	  write(11,*) "Value   Contacts  (sp1site-sp2site,...)"
	  do n=1,2**(2*nsp2sites(sp2))-1
	    write(strtmp,"(i4,2x)") n
	    i = 7
	    do site=1,nsp2sites(sp2)
	      !if (btest(n,site-1)) write(0,"(i2,'-',i2)") sp1sites(1), sp2sites(sp2,site,1)
	      !if (btest(n,nsp2sites(sp2)+site-1)) write(0,"(i2,'-',i2)") sp1sites(2), sp2sites(sp2,site,1)
	      if (btest(n,site-1)) then
		strtmp = strtmp(1:i)//CHAR(48+(sp1sites(1)/10))//CHAR(48+MOD(sp1sites(1),10))//"-" &
		& //CHAR(48+(sp2sites(sp2,site,1)/10))//CHAR(48+MOD(sp2sites(sp2,site,1),10))//" "
		i = i + 6
	      end if
	      if (btest(n,nsp2sites(sp2)+site-1)) then
		strtmp = strtmp(1:i)//CHAR(48+(sp1sites(2)/10))//CHAR(48+MOD(sp1sites(2),10))//"-" &
		& //CHAR(48+(sp2sites(sp2,site,1)/10))//CHAR(48+MOD(sp2sites(sp2,site,1),10))//" "
		i = i + 6
	      end if
	    end do
	    write(11,*) strtmp(1:i)
	  end do
	  write(11,*) ""
	end do

	! Output pattern data
	write(11,*) "Number of distinct contact patterns found:",npatterns
	write(11,*) ""
	write(11,"(a)") "  ID       Freq       %     Cumul%  Nm  Nc  Pattern"
850	FORMAT (i4,2x,f12.8,2x,f6.2,2x,f6.2,2x,i2,2x,i2,2x,20(i2,'/',i6,','))
	avgatom = 0.0
	avgmol = 0.0
	total = 0.0
	cp = 0.0
	do n=1,npatterns
	  cp = cp + 100.0*pfreq(n)/(s_nmols(sp1)*nframes)
	  write(11,850) n, pfreq(n)*1.0/nframes, 100.0*pfreq(n)/(s_nmols(sp1)*nframes), cp, nmolcontacts(n), natomcontacts(n), (pattern(n,i,1),pattern(n,i,3),i=1,nmolcontacts(n)+1)
	  total = total + pfreq(n)*1.0/nframes
	  avgmol = avgmol + nmolcontacts(n) * pfreq(n)*1.0/nframes
	  avgatom = avgatom + natomcontacts(n) * pfreq(n)*1.0/nframes
	end do
851	FORMAT ('SNTY',2x,f12.4,7x,2(f4.2,1x))
	write(11,851) total, avgmol / s_nmols(sp1), avgatom / s_nmols(sp1)

	write(11,*) ""
	write(11,"(a,i5,a)") "Total contacts found for defined sites over ", nframes," frames"
	do sp2=1,nspecies
	  do site=1,nsp2sites(sp2)
	    write(11,"(2i4,2i8,', per frame = ',f12.5)") sp2, site, nsp2atomcontacts(sp2,site,1), nsp2atomcontacts(sp2,site,2),sum(nsp2atomcontacts(sp2,site,1:2))*1.0/nframes
	  end do
	end do

	close(11)

	! Write out 3d data if requested
	if (calc3d) then
	  do n=1,npatterns
	    write(0,*) "Writing 3D data for pattern ",n," of ",npatterns
	    ! Determine number of atoms to write
	    natoms = s_natoms(sp1)
	    do m=1,nmolcontacts(n)
	      natoms = natoms + s_natoms(pattern(n,m,1))
	    end do
	    ! Normalise coordinates
	    if (.not.writeone) geom(n,:,:) = geom(n,:,:) / (pfreq(n)*1.0)
	    ! Construct destfile name
	    tth = n / 10000; i = n - tth*10000
	    th = i / 1000; i = i - th*1000
	    hun = i / 100; i = i - hun*100
	    ten = i / 10; i = i - ten*10
	    units = i
	    open(unit=12,file=basename(1:baselen)//"bident4_"//CHAR(48+sp1)//"_p"//CHAR(48+tth)//CHAR(48+th)//CHAR(48+hun)//CHAR(48+ten)//CHAR(48+units)//".xyz",form='formatted',status='replace')
	    write(12,*) natoms
	    write(12,*) "Pattern ",n
	    aoff1 = s_natoms(sp1)
	    do i=1,s_natoms(sp1)
	      write(12,"(a8,2x,3f12.6)") atmname(s_start(sp1)-1+i), geom(n,i,1), geom(n,i,2), geom(n,i,3)
	    end do
	    do m=1,nmolcontacts(n)
	      sp2 = pattern(n,m,1)
	      aoff2 = s_start(sp2) - 1
	      do i=1,s_natoms(sp2)
	        write(12,"(a8,2x,3f12.6)") atmname(aoff2+i), geom(n,aoff1+i,1), geom(n,aoff1+i,2), geom(n,aoff1+i,3)
	      end do
	      aoff1 = aoff1 + s_natoms(sp2)
	    end do
	    close(12)
	    open(unit=12,file=basename(1:baselen)//"bident4_"//CHAR(48+sp1)//"_p"//CHAR(48+tth)//CHAR(48+th)//CHAR(48+hun)//CHAR(48+ten)//CHAR(48+units)//".pdens",form='formatted',status='replace')
	    write(12,*) 2*grid+1, 2*grid+1, 2*grid+1
	    write(12,"(9f8.4)") delta,0.0,0.0,0.0,delta,0.0,0.0,0.0,delta
	    write(12,"(3f10.4)") -grid*delta,-grid*delta,-grid*delta
	    write(12,*) "zyx"
	    do i=-grid,grid
	      do j=-grid,grid
	        do k=-grid,grid
	          write(12,"(f12.8)") pdens(n,i,j,k) / pfreq(n)
	        end do
	      end do
	    end do
	    close(12)

	  end do
	end if

	write(0,*) "Finished."
999	close(14)
	end program bident4
