!	** cluster **
!	Calculate cluster analysis between COMS of molecules

	program cluster
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	logical :: altheader = .FALSE.
	integer :: n,sp1,m1,m2,baselen,nframes,success,nargs,framestodo = -1,frameskip = 0,framesdone,compairs(10,2)
	integer, allocatable :: marked(:)
	integer :: iargc, nmarked, oldnmarked, newsize
	real*8 :: c1x,c1y,c1z,c2x,c2y,c2z,tx,ty,tz,maxdist
	real*8, allocatable :: dist(:,:), clusters(:)

	compairs = 0
	nargs = iargc()
	if (nargs.lt.4) stop "Usage : cluster <DLP HISTORYfile> <DLP OUTPUTfile> <sp> <maxcomdist> [-header hisfile] [-frames n] [-discard n] [-compair sp i j]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp)
        read(temp,"(i10)") sp1
	call getarg(4,temp)
        read(temp,"(f10.4)") maxdist
	
	n = 4
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A,I4)") "Alternative header file supplied."
	      altheader = .TRUE.
            case ("-frames")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
              write(0,"(A,I4)") "Frames to process: ",framestodo
            case ("-discard")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") frameskip
              write(0,"(A,I4)") "Frames to discard at start: ",frameskip
            case ("-compair")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") sp1
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(sp1,1)
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(sp1,2)
              write(0,"(A,3I4)") "Using COMpair for species ",sp1, compairs(sp1,:)
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
	allocate(dist(s_nmols(sp1),s_nmols(sp1)))
	allocate(marked(s_nmols(sp1)))
	allocate(clusters(s_nmols(sp1)))
	dist = 0.0
	marked = 0
	clusters = 0.0
	nmarked = 0

	! Set up the vars...
100	nframes=0
	framesdone = 0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.EQ.-1) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes, framesdone
	if (nframes.le.frameskip) goto 101

	framesdone = framesdone + 1

	call calc_com

	! If compairs were specified, use that instead of COM
	if (compairs(sp1,1).ne.0) then
	  do m1=1,s_nmols(sp1)
	    c1x = xpos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,1)-1)
	    c1y = ypos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,1)-1)
	    c1z = zpos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,1)-1)
	    c2x = xpos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,2)-1)
	    c2y = ypos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,2)-1)
	    c2z = zpos(s_start(sp1)+(m1-1)*s_natoms(sp1)+compairs(sp1,2)-1)
	    call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	    comx(sp1,m1) = (c1x+tx)*0.5
	    comy(sp1,m1) = (c1y+ty)*0.5
	    comz(sp1,m1) = (c1z+tz)*0.5
	  end do
	end if

	! Calculate distance array
	dist = 0.0
	do m1=1,s_nmols(sp1)-1     ! Loop over all molecules of species 1...
	  c1x=comx(sp1,m1)
	  c1y=comy(sp1,m1)
	  c1z=comz(sp1,m1)
	  ! Now loop over all molecules of second species....
	  do m2=m1+1,s_nmols(sp1)
	    ! Grab the centre-of-mass coordinates...
	    c2x=comx(sp1,m2)
	    c2y=comy(sp1,m2)
	    c2z=comz(sp1,m2)
	    ! Get the shortest (MIM) distance between the atom pair...
	    call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	    dist(m1,m2)=sqrt( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )
	    dist(m2,m1) = dist(m1,m2)
	  end do
	end do

	! Calculate cluster sizes
	marked = 0
	nmarked = 0
	do m1=1,s_nmols(sp1)
	  if (marked(m1).eq.1) cycle
	  oldnmarked = nmarked
	  ! From this molecule, mark it and all neighbours within maxdist
	  call mark(m1, s_nmols(sp1), nmarked, marked, maxdist, dist)
	  newsize = nmarked - oldnmarked
	  !write(0,*) "Cluster size is ", newsize
	  clusters(newsize) = clusters(newsize) + 1
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
	write(9,"(a,i2,a,f10.4)") "# Species ",sp1," clusters with COM distance < ", maxdist
	if (compairs(sp1,1).eq.0) then
	  write(9,"(a,i2,a)") "# Species ",sp1," calculated using all atoms for COM."
	else
	  write(9,"(a,i2,a,i3,i3)") "# Species ",sp1," calculated using two atoms for COM: ", compairs(sp1,1), compairs(sp1,2)
	end if
	! Normalise the RDFs with respect to the number of frames.
	clusters = clusters / framesdone
	do n=1,s_nmols(sp1)
	  ! write(9,"(F6.3,3x,F12.8)") (n*binwidth)-binwidth/2.0,nrdf(a1,n)
	  write(9,"(i6,2(3x,F12.8))") n, clusters(n)
	end do
	close(9)
	write(0,*) "Finished."
999	close(10)
	close(13)

	! Deallocate arrays
	deallocate(dist)
	deallocate(marked)
	deallocate(clusters)

	end program cluster

	recursive subroutine mark(mol, nmols, nmarked, marked, maxdist, dist)
	implicit none
	integer, intent(inout) :: nmarked, marked(nmols)
	integer, intent(in) :: mol, nmols
	integer :: m1
	real*8, intent(in) :: maxdist, dist(nmols,nmols)
	! Mark the 'mol' specified, if it has not been marked already
	if (marked(mol).eq.0) then
	  marked(mol) = 1
	  nmarked = nmarked + 1
	end if
	! Find close neighbours, and mark those as well
	do m1=1,nmols
	  if ((m1.eq.mol).or.(dist(mol,m1).gt.maxdist)) cycle
	  ! Within distance range, so mark
	  if (marked(m1).eq.0) then
	    marked(m1) = 1
	    nmarked = nmarked + 1
	    call mark(m1, nmols, nmarked, marked, maxdist, dist)
	  end if
	end do
	end subroutine mark
