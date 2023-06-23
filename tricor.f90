	program tricor
	use dlprw; use utility; use parse
	implicit none
	! File names
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	! Calculation Control
	logical :: altheader = .FALSE.
	integer :: framestodo = -1,frameskip = 0,framesdone
	! Correlation matrix data
	integer :: nbins_ab, nbins_cd, bin_ab, bin_cd
	real*8 :: binwidth
	real*8, allocatable :: cormat_cd(:,:)
	! Working variables
	real*8, parameter :: pi = 3.14159265358979d0
	integer :: iargc
	integer :: numadded,n,m,i,a,b,c,d,sp1,sp2,sp3,m1,m2,m3,aoff1,aoff2,aoff3,baselen,bin,nframes,success,o,nargs, status
	real*8 :: dist,ra(3), rb(3), rc(3), rd(3), tmp(3),bdens,integral, cellvol, norm, rmax_ab, rmax_bc, rmax_cd = 10.0
	real*8 :: rab, rbc, rcd

	binwidth=0.05   ! In Angstroms
	nargs = iargc()
	if (nargs.LT.11) stop "Usage : tricor <HISTORYfile> <OUTPUTfile> <sp1> <sp2> <sp3> <a> <b> <rmax_ab> <c> <rmax_bc> <d> [-bin width] [-header hisfile] [-frames n] [-discard n]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	sp1 = getargi(3)
	sp2 = getargi(4)
	sp3 = getargi(5)
	a = getargi(6)
	b = getargi(7)
	rmax_ab = getargi(8)
	c = getargi(9)
	rmax_bc = getargi(10)
	d = getargi(11)
	
	n = 11
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-bin")
              n = n + 1; binwidth = getargr(n)
              write(0,"(A,f6.2)") "Binwidth set to ",binwidth
            case ("-discard")
              n = n + 1; frameskip = getargi(n)
              write(0,"(A,I4)") "Frames to discard at start: ",frameskip
            case ("-frames")
              n = n + 1; framestodo = getargi(n)
              write(0,"(A,I4)") "Frames to process: ",framestodo
            case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A,I4)") "Alternative header file supplied."
	      altheader = .TRUE.
	    case default
	      write(0,"(a,a)") "Unrecognised command line option:",temp
	      stop
	  end select
	end do

	! Open and check the files...
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798

        ! Now, read in the history header so that we have cell()
	! If this fails then we may have a restarted trajectory. Continue, but only if
	! a header can be read in from the specified alternative history file..
	call openhis(hisfile,10)
        if (readheader().EQ.-1) then
	  if (altheader) then
	    write(0,*) "Restarted trajectory:"
	    close(dlpun_his)
	    call openhis(altheaderfile,10)
	    if (readheader().EQ.-1) goto 797
	    close(dlpun_his)
	    call openhis(hisfile,10)
	  else
	    goto 797
	  end if
	end if

	nbins_ab = rmax_ab / binwidth + 1
	nbins_cd = rmax_cd / binwidth + 1
	
	! Allocate arrays
        allocate(cormat_cd(nbins_ab, nbins_cd),stat=status); if (status.GT.0) stop "Allocation error for cormat_cd()"

	! Initialise the arrays...
	cormat_cd = 0.0

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
100	nframes=0
	framesdone = 0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes, framesdone
	if (nframes.le.frameskip) goto 101

	framesdone = framesdone + 1

	! Outer loop over molecules of species 1 (containing atoms 'a' and 'd')
	aoff1 = s_start(sp1) - 1 - s_natoms(sp1)
	do m1=1,s_nmols(sp1)

	  aoff1 = aoff1 + s_natoms(sp1)

	  ! Grab coordinates of atoms 'a' and 'd' (the latter MIM'd with the first)
	  ra(1) = xpos(aoff1+a)
	  ra(2) = ypos(aoff1+a)
	  ra(3) = zpos(aoff1+a)
	  tmp(1) = xpos(aoff1+d)
	  tmp(2) = ypos(aoff1+d)
	  tmp(3) = zpos(aoff1+d)
	  call pbc(tmp(1), tmp(2), tmp(3), ra(1), ra(2), ra(3), rd(1), rd(2), rd(3))

	  ! Middle loop over molecules of species 2 (containing atom 'b')
	  aoff2 = s_start(sp2) - 1 - s_natoms(sp2)
	  do m2=1,s_nmols(sp2)

	    aoff2 = aoff2 + s_natoms(sp2)

	    ! Grab coordinates of 'b', mim'd w.r.t. 'a'
	    tmp(1) = xpos(aoff2+b)
	    tmp(2) = ypos(aoff2+b)
	    tmp(3) = zpos(aoff2+b)
	    call pbc(tmp(1), tmp(2), tmp(3), ra(1), ra(2), ra(3), rb(1), rb(2), rb(3))

	    ! Check distance a-b
	    tmp = ra - rb
	    rab = sqrt(tmp(1)*tmp(1) + tmp(2)*tmp(2) + tmp(3)*tmp(3))
	    if (rab.gt.rmax_ab) cycle

	    ! Distance a-b is OK, so find any molecules of sp3 that are in range...
	    bin_ab =  rab / binwidth
	    aoff3 = s_start(sp3) - 1 - s_natoms(sp3)
	    do m3=1,s_nmols(sp3)

	      aoff3 = aoff3 + s_natoms(sp3)

	      ! Grab coordinates of 'c', mim'd w.r.t. 'b'
	      tmp(1) = xpos(aoff3+c)
	      tmp(2) = ypos(aoff3+c)
	      tmp(3) = zpos(aoff3+c)
	      call pbc(tmp(1), tmp(2), tmp(3), rb(1), rb(2), rb(3), rc(1), rc(2), rc(3))

	      ! Check distance b-c
	      tmp = rb - rc
	      rbc = sqrt(tmp(1)*tmp(1) + tmp(2)*tmp(2) + tmp(3)*tmp(3))
	      if (rbc.gt.rmax_bc) cycle

	      ! AB and BC are both in range, so get mim'd 'c' w.r.t. 'd' (to be on the safe side) and add to cormat_cd
	      call pbc(rc(1), rc(2), rc(3), rd(1), rd(2), rd(3), rc(1), rc(2), rc(3))
	      tmp = rc - rd
	      rcd = sqrt(tmp(1)*tmp(1) + tmp(2)*tmp(2) + tmp(3)*tmp(3))
	      if (rcd.gt.rmax_cd) then
		write(0,*) "Missed one at ", rcd
		cycle
	      end if
	      bin_cd = rcd / binwidth
	      
	      cormat_cd(bin_ab, bin_cd) = cormat_cd(bin_ab, bin_cd) + 1.0

	    end do ! m3 of sp3

	  end do ! m2 of sp2

	end do ! m1 of sp1

	if (framesdone.eq.framestodo) goto 801
	! Next frame
	goto 101

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

797	write(0,*) "No header found in history file. If a restarted trajectory, use '-header'"
	goto 999
798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended prematurely!"
	goto 801
800	write(0,*) "End of unformatted HISTORY file found."
801	write(0,*) ""

	! Normalise and write out data
	cormat_cd = cormat_cd / nframes
	resfile=outputFileName(hisfile, "corcd", stringNMMOO(sp1,a,b)//"_"//stringNMM(sp2,c)//"_"//stringNMM(sp3,d)//".corcd")
	open(unit=9,file=resfile,form="formatted")
	do n=1,nbins_ab
	  do m=1,nbins_cd
	    write(9,*) cormat_cd(n,m)
	  end do
	  write(9,*) ""
	end do
	close(9)
	
999	write(0,*) "Finished."
	end program tricor

