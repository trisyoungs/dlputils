!	** rdf_aa_inter **
!	Calculate RDFs between specific atoms of two species

	program rdf_ss
	use dlprw; use utility
	implicit none
	integer, parameter :: maxpairs = 50
	real*8, parameter :: pi = 3.14159265358979d0
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	integer :: status,  nargs, success, baselen, npairs
	integer :: nbins,aoff1,aoff2,n,s1,s2,bin,nframes,sp1,sp2,atom1(maxpairs),atom2(maxpairs)
	integer :: iargc,p,framestodo,framestodiscard=0
	logical :: nonorm = .FALSE., includeintra = .FALSE.
	real*8 :: r1(3), r2(3), r2min(3), r12(3), dist, binwidth, norm, integral, dumpdist = -1.0, numadded
	real*8, allocatable :: hist(:,:), rdf(:,:), sumhist(:,:)

	binwidth=0.01   ! In Angstroms
	framestodo=-1
	sp1 = 1
	sp2 = 2
	npairs = 0

	nargs = iargc()
	if (nargs.lt.5) stop "Usage : rdfaainter <HISTORYfile> <OUTPUTfile> [-sp1 n] [-sp2 n] -pair a1 a2 [-pair a1 a2 [...] ] [-frames n] [-nonorm] [-dump <dist>] [-includeintra] [-bin width]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	n = 2
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-bin") 
	      n = n + 1; call getarg(n,temp); read(temp,"(f20.14)") binwidth
	      write(0,"((A),f9.6)") "Bin width set to ", binwidth
	    case ("-sp1") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp1
	    case ("-sp2") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp2
	    case ("-frames") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I5)") framestodo
	    case ("-discard") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I5)") framestodiscard
	    case ("-pair")
	      npairs = npairs + 1
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") atom1(npairs)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") atom2(npairs)
            case ("-nonorm")
              write(0,"(A)") "RDFs will not be normalised."
	      nonorm = .TRUE.
            case ("-dump")
	      n = n + 1; call getarg(n,temp); read(temp,"(f12.5)") dumpdist
              write(0,"(A)") "Dump contacts with distance less than", dumpdist
	    case ("-includeintra")
	      write(0,"(A)") "Intramolecular contacts will be included (if sp1 == sp2)."
	      includeintra = .TRUE.
	    case default
	      write(0,*) "Unrecognised option: ", temp
	      stop
	  end select
	end do
	     
	if (npairs.eq.0) stop "No atom pairs specified!"

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	nbins = maxval(cell) / binwidth
	write(0,"(A,i1,a,i1)") "'a1 / a2 taken from species ",sp1," / ",sp2
	write(0,"(A,F6.3,A)") "Using binwidth of ",binwidth," Angstroms"
	write(0,"(A,I5,A)") "There will be ",nbins," histogram bins."
	if (framestodo.gt.0) write(0,"(a,i6)") "Number of frames to use in average = ",framestodo
	write(0,"(a,20(i3,'/',i3,','))") "Pairs to calculate are : ",(atom1(n),atom2(n),n=1,npairs)
	
        allocate(rdf(npairs,nbins),stat=status); if (status.GT.0) stop "Allocation error for rdf()"
	allocate(hist(npairs,nbins),stat=status); if (status.GT.0) stop "Allocation error for nrdf()"
	allocate(sumhist(npairs,nbins),stat=status); if (status.GT.0) stop "Allocation error for nrdf()"

	! Initialise the arrays...
	rdf = 0.0
	sumhist = 0

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
	numadded = 0.0
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 120  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	if (framestodiscard.gt.0) then
	  framestodiscard = framestodiscard -1
	  goto 101
	end if
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	
	hist = 0.0

	aoff1 = s_start(sp1)
	do s1 = 1,s_nmols(sp1)     ! Loop over all molecules in species1

	  aoff2 = s_start(sp2) 
	  do s2 = 1,s_nmols(sp2)     ! Loop over all molecules of species 2

	    ! Check for intramolecular exclusion
	    if ((sp1.eq.sp2).and.(s1.eq.s2).and.(.not.includeintra)) then
	      aoff2 = aoff2 + s_natoms(sp2)
	      cycle
	    end if

	    ! Loop over all pairs of atoms specified
	    do p = 1,npairs

	      ! Grab coordinates of what will be the central atom
	      r1(1) = xpos(aoff1+atom1(p)-1)
	      r1(2) = ypos(aoff1+atom1(p)-1)
	      r1(3) = zpos(aoff1+atom1(p)-1)
  
	      r2(1) = xpos(aoff2+atom2(p)-1)
	      r2(2) = ypos(aoff2+atom2(p)-1)
	      r2(3) = zpos(aoff2+atom2(p)-1)

	      call pbc(r2(1),r2(2),r2(3),r1(1),r1(2),r1(3),r2min(1),r2min(2),r2min(3))
	      r12 = r2min - r1
	      dist = sqrt( r12(1)*r12(1) + r12(2)*r12(2) + r12(3)*r12(3) )
  
	      if (dist.lt.dumpdist) write(66,*) aoff1+atom1(p)-1, aoff2+atom2(p)-1, dist
	      bin = int(dist * (1.0/binwidth))+1
	      hist(p,bin) = hist(p,bin)+1

	    end do
  
	    numadded = numadded+1.0
	    aoff2 = aoff2 + s_natoms(sp2)
	  end do

	  aoff1 = aoff1 + s_natoms(sp1)

	end do   ! End main loop over all atoms of species1.

	if (nframes.EQ.1) then
	  write(0,"(A,I2,A,I4,A)") "PRDF of atom in ",sp1," : averaged over ",s_nmols(sp1)," molecules."
	end if

	! Sum histogram into final RDF array
	sumhist = sumhist + hist

	! Next frame
	if ((framestodo.gt.0).and.(framestodo.eq.nframes)) goto 120
	goto 101
	! Work on the results now to get the proper RDF
120	write(0,*) "Finished."
	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended prematurely!"
	goto 801
800	write(0,*) "End of unformatted HISTORY file found."

801	write(0,*) ""
	write(0,"(a,i6)") "Total frames used = ",nframes
	write(0,"(a,e12.5,a,e12.5)") "Total numadded = ",numadded,", per frame = ",numadded/nframes
	write(0,"(a,e12.5)") "Total numadded per pair ",numadded / npairs

	! Normalise data
	do n=1,nbins
	  norm = (4.0 * pi / 3.0) * ((n*binwidth)**3 - ((n-1)*binwidth)**3) / (cell(1)*cell(5)*cell(9))
	  do p=1,npairs
	    if (nonorm) then
	      rdf(p,n) = sumhist(p,n) / nframes
	    else
	      rdf(p,n) = sumhist(p,n) / norm / s_nmols(sp1) / s_nmols(sp2) / nframes
	    end if
	  end do
	end do

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).EQ.".") THEN
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) THEN
	  basename="rdfresults."
	  baselen=11
	ELSE
	  basename=hisfile(1:baselen)
	endif

	do p=1,npairs
	  resfile=basename(1:baselen)//"aardf"//CHAR(48+sp1)//"_"//CHAR(48+(atom1(p)/10))//CHAR(48+MOD(atom1(p),10))// &
		& "_"//CHAR(48+sp2)//"_"//CHAR(48+(atom2(p)/10))//CHAR(48+MOD(atom2(p),10))
	  OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	  integral = 0.0
	  do n=1,nbins
	    integral = integral + sumhist(p,n) / nframes / s_nmols(sp1)
	    ! write(0,"(i4,2x,e12.6,2x,e12.6,2x,e12.6)") n,rdf(p,n),sumhist(p,n)
! 	    write(9,"(f6.3,3x,f12.8,4x,e12.6,2x,e12.6)") (n*binwidth)-binwidth/2.0,rdf(p,n),sumhist(p,n),integral
	    write(9,"(f10.4,3x,f12.8,4x,e12.6,2x,e12.6)") (n*binwidth)-binwidth/2.0,rdf(p,n),integral,sumhist(p,n)
	  end do
	  close(9)
	end do

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end

