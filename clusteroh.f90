!	** clusteroh **
!	Calculate cluster analysis between COMS of molecules

	program clusteroh
	use parse; use dlprw; use utility; use IList
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	logical :: altheader = .FALSE., individual = .false., notSelf = .true.
	integer :: n,sp1,sp,m1,m2,baselen,nframes,success,nargs,framestodo = -1,frameskip = 0,framesdone
	type(IntegerList) :: otherSp, oxygenAtoms(MAXSP), hydrogenAtoms(MAXSP)
	integer, allocatable :: marked(:)
	integer :: iargc, newsize, oldsize, maxpath = -1
	real*8 :: maxdist, sumtotal
	real*8, allocatable :: clusters(:)

	nargs = iargc()
	if (nargs.lt.5) stop "Usage : clusteroh <HISTORYfile> <OUTPUTfile> <sp> <othersp> <maxdist> [-oh sp Oatoms Hatoms] [-header hisfile] [-frames n] [-discard n] [-individual] [-maxpath n]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
        sp1 = getargi(3)
	call getarg(4,temp); if (.not.parseIntegerList(temp, otherSp)) stop "Failed to parse other species list." 
        maxdist = getargr(5)
	
	n = 5
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-oh")
	      n = n + 1; sp = getargi(n)
	      n = n + 1; call getarg(n,temp); if (.not.parseIntegerList(temp, oxygenAtoms(sp))) stop "Failed to parse oxygen atom list."
	      n = n + 1; call getarg(n,temp); if (.not.parseIntegerList(temp, hydrogenAtoms(sp))) stop "Failed to parse hydrogen atom list."
	      write(0,"(i2,a,i2,a,20i4)") oxygenAtoms(sp)%n, " oxygen atom sites added as for species ", sp, ": ", oxygenAtoms(sp)%items(1:oxygenAtoms(sp)%n)
	      write(0,"(i2,a,i2,a,20i4)") hydrogenAtoms(sp)%n, " hydrogen atom sites added as for species ", sp, ": ", hydrogenAtoms(sp)%items(1:hydrogenAtoms(sp)%n)
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
            case ("-individual")
	      individual = .true.
              write(0,"(A)") "Individual mode - marked list will be reset for each molecule of target species."
	    case ("-maxpath")
	      n = n + 1; maxpath = getargi(n)
	      write(0,"(a,i3)") "Maximum path for selection from each molecule of target species set to ", maxpath
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

	! Allocate arrays
	allocate(marked(s_totalMols))
	allocate(clusters(s_totalMols))
	clusters = 0.0
	notSelf = .not.listContains(otherSp,sp1)

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

	! Loop over each molecule of sp1, marking from each
	! Zero marked array each time if -individual was specified
	marked = 0

	do m1=1,s_nmols(sp1)
	  ! Cycle if this molecule is already marked
	  if (marked(s_molstart(sp1)+m1-1).eq.1) cycle

	  ! From this molecule, mark it and all neighbours within maxdist
	  oldsize = sum(marked)
	  call mark(sp1, m1, marked, maxdist, otherSp, oxygenAtoms, hydrogenAtoms, 0, maxPath)
	  newsize = sum(marked) - oldsize
	  !write(0,*) "Cluster size is ", newsize
	  clusters(newsize) = clusters(newsize) + 1

	  ! Zero array?
	  if (individual) marked = 0
	end do

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

	! Write output
	resfile=basename(1:baselen)//"cluster"//CHAR(48+sp1)
	open(unit=9,file=resfile,form="formatted")
	!write(9,"(a,i2,a,f10.4)") "# Species ",sp1," clusters with COM distance < ", maxdist

	! Normalise w.r.t. number of frames
	clusters = clusters / framesdone

	! Write data
	write(9,"(a6,5(3x,a12))") "# Size", "NClusters", "NMolecules", "%Clusters", "%Molecules", "Sum(S*N)"
	sumtotal = 0.0
	if (notSelf) then
	  do n=1,s_nmols(sp1)
	    sumtotal = sumtotal + (n-1)*clusters(n)
	    write(9,"(i6,5(3x,F12.8))") n-1, clusters(n), n*clusters(n), 100.0*real(clusters(n))/sum(clusters), 100.0*n*clusters(n)/s_nmols(sp1), sumtotal
	  end do
	else
	  do n=1,s_nmols(sp1)
	    sumtotal = sumtotal + n*clusters(n)
	    write(9,"(i6,5(3x,F12.8))") n, clusters(n), n*clusters(n), 100.0*real(clusters(n))/sum(clusters), 100.0*n*clusters(n)/s_nmols(sp1), sumtotal
	  end do
	end if
	close(9)
	write(0,*) "Finished."
999	close(10)
	close(13)

	! Deallocate arrays
	deallocate(marked)
	deallocate(clusters)

	end program clusteroh

	recursive subroutine mark(currentSp, currentMol, marked, maxdist, otherSp, oxygenAtoms, hydrogenAtoms, pathSize, maxPath)
	use dlprw; use utility; use IList
	implicit none
	integer, intent(inout) :: marked(s_totalMols)
	integer, intent(in) :: currentSp, currentMol, pathSize, maxPath
	integer :: currentAtomOffset, currentMolOffset
	integer :: moff, m, sp, mol, aoff, i, j, newPathSize
	type(IntegerList), intent(in) :: otherSp, oxygenAtoms(MAXSP), hydrogenAtoms(MAXSP)
	real*8, intent(in) :: maxdist
	real*8 :: O(3), H(3), v(3), dist

	! Calculate atom and molecule offsets for the current molecule / species
	currentAtomOffset = (s_start(currentSp)-1) + (currentMol-1)*s_natoms(currentSp)
	currentMoloffset = s_molstart(currentSp)-1

	! Mark the 'mol' specified, if it has not been marked already
	if (marked(currentMol+currentMolOffset).eq.0) then
	  marked(currentMol+currentMolOffset) = 1
	end if

	! Increase pathsize
	newPathSize = pathSize + 1

	! Find close neighbours (over active species), and mark those as well
	aoff = 0
	do sp=1,nspecies
	  ! Loop if this species is not in the otherSp list
	  if (.not.listContains(otherSp,sp)) cycle

	  ! Loop over all molecules of this species, checking for H...O and O...H contacts with the current molecule
	  aoff = s_start(sp)-1
	  moff = s_molstart(sp)-1
	  do m=1,s_nmols(sp)

	    ! Skip if this molecule/species are the same as the current, or if it is already selected
	    if ((sp.eq.currentSp).and.(m.eq.currentMol)) then
	      aoff = aoff + s_natoms(sp)
	      cycle
	    end if
	    if (marked(m+moff).eq.1) then
	      aoff = aoff + s_natoms(sp)
	      cycle
	    end if

	    ! Check O (on currentMol) to H (on our molecule)
	    do i=1,oxygenAtoms(currentSp)%n

	      ! Grab coordinates of the oxygen
	      O(1) = xpos(currentAtomOffset + oxygenAtoms(currentSp)%items(i))
	      O(2) = ypos(currentAtomOffset + oxygenAtoms(currentSp)%items(i))
	      O(3) = zpos(currentAtomOffset + oxygenAtoms(currentSp)%items(i))

	      ! Loop over H sites on our molecule
	      do j=1,hydrogenAtoms(sp)%n

	        ! Grab coordinates of the hydrogen
		v(1) = xpos(aoff + hydrogenAtoms(sp)%items(j))
		v(2) = ypos(aoff + hydrogenAtoms(sp)%items(j))
		v(3) = zpos(aoff + hydrogenAtoms(sp)%items(j))
		
		! Get mim coordinates of H w.r.t. O
		call pbc(v(1),v(2),v(3),O(1),O(2),O(3),H(1),H(2),H(3))

		! Get and check distance
		v = H - O
		dist = dsqrt(sum(v*v))
		if (dist.gt.maxdist) cycle

		! Within range, so mark it (and possibly those nearby...)
		if ((maxPath.gt.0).and.(newPathSize.eq.maxPath)) then
		  marked(m+moff) = 1
		else
		  call mark(sp, m, marked, maxdist, otherSp, oxygenAtoms, hydrogenAtoms, newPathSize, maxPath)
		end if

	      end do !H

	    end do !O

	    ! Check H (on currentMol) to O (on our molecule)
	    do i=1,hydrogenAtoms(currentSp)%n

	      ! Grab coordinates of the hydrogen
	      H(1) = xpos(currentAtomOffset + hydrogenAtoms(currentSp)%items(i))
	      H(2) = ypos(currentAtomOffset + hydrogenAtoms(currentSp)%items(i))
	      H(3) = zpos(currentAtomOffset + hydrogenAtoms(currentSp)%items(i))

	      ! Loop over O sites on our molecule
	      do j=1,oxygenAtoms(sp)%n

	        ! Grab coordinates of the oxygen
		v(1) = xpos(aoff + oxygenAtoms(sp)%items(j))
		v(2) = ypos(aoff + oxygenAtoms(sp)%items(j))
		v(3) = zpos(aoff + oxygenAtoms(sp)%items(j))
		
		! Get mim coordinates of O w.r.t. H
		call pbc(v(1),v(2),v(3),H(1),H(2),H(3),O(1),O(2),O(3))

		! Get and check distance
		v = H - O
		dist = dsqrt(sum(v*v))
		if (dist.gt.maxdist) cycle

		! Within range, so mark it (and possibly those nearby...)
		if ((maxPath.gt.0).and.(newPathSize.eq.maxPath)) then
		  marked(m+moff) = 1
		else
		  call mark(sp, m, marked, maxdist, otherSp, oxygenAtoms, hydrogenAtoms, newPathSize, maxPath)
		end if

	      end do !H

	    end do !O

	    ! Increase atom offset and continue
	    aoff = aoff + s_natoms(sp)

	  end do

	end do

	end subroutine mark
