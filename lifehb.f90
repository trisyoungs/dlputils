!	** lifehb **
!	Determine H-Bond lifetime histogram for interactions specified

	program lifehb
	use dlprw; use utility
	implicit none
	integer, parameter :: MAXXH = 20, MAXY = 20, MAXBINS = 10000
	character*80 :: hisfile,dlpoutfile,altheaderfile
	character*80 :: temp
	integer :: n,nframes,success,nargs,m,i,iargc,m1,m2,j
	integer :: framestodo = -1, framestodiscard = 0, t(9), sp1, sp2, nxh = 0, ny = 0, framesdone = 0
	integer :: xh(MAXXH,2), y(MAXY), ncontacts, maxlife = 0
	integer, allocatable :: hbflag(:,:,:,:)
	real*8 :: dxy, rx(3), rh(3), vhx(3), dhx, ry(3), vhy(3), dhy, angle, dp, maxdist, minang, histo(MAXBINS)
	logical :: altheader = .FALSE.

	nargs = iargc()
	if (nargs.LT.7) then
	  write(0,"(a)") "Usage : lifehb <HISTORYfile> <OUTPUTfile> <sp1 (XH)> <sp2 (Y)> <rXY max> <aXHY min>"
	  write(0,"(10x,a)") "[-xh i j ...]     Specify donor site on species 1"
	  write(0,"(10x,a)") "[-y k ...]        Specify acceptor site on species 2"
	  write(0,"(10x,a)") "[-frames n]       Set number of frames to calculate (default = all)"
	  write(0,"(10x,a)") "[-discard n]      Set number of frames to discard at start (default = 0)"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(i4)") sp1
	call getarg(4,temp); read(temp,"(i4)") sp2
	call getarg(5,temp); read(temp,"(f12.4)") maxdist
	call getarg(6,temp); read(temp,"(f12.4)") minang

	n = 6
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
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodiscard
              write(0,"(A,I4)") "Frames to discard at start: ",framestodiscard
            case ("-xh")
	      nxh = nxh + 1
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") xh(nxh,1)
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") xh(nxh,2)
            case ("-y")
	      ny = ny + 1
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") y(ny)
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
	    write(6,*) "Restarted trajectory:"
	    close(dlpun_his)
	    call openhis(altheaderfile,10)
	    if (readheader().EQ.-1) goto 799
	    close(dlpun_his)
	    call openhis(hisfile,10)
	  else
	    goto 799
	  end if
	end if

	write(6,*) "Species 1 (XH) ", sp1, s_name(sp1)
	write(6,*) "Species 2 (Y)  ", sp2, s_name(sp2)
	write(6,*) "Maximum X...Y distance :",maxdist
	write(6,*) "Minimum X-H...Y angle  :",minang
	if (nxh.eq.0) stop "Error: No XH sites defined on species 1"
	write(6,*) "Number of XH sites defined on species 1: ", nxh
	write(6,"(12x,i3,' (',a8,') - ',i3,' (',a8,')')") (xh(n,1), s_atom(sp1,xh(n,1)), xh(n,2), s_atom(sp1,xh(n,2)), n=1,nxh)
	if (nxh.eq.0) stop "Error: No Y sites defined on species 2"
	write(6,*) "Number of Y sites defined on species 2: ", ny
	write(6,"(12x,i3,' (',a8,')')") (y(n), s_atom(sp2,y(n)), n=1,ny)

	! Allocate lifetime flag array
	allocate(hbflag(s_nmols(sp1),nxh,s_nmols(sp2),ny))
	hbflag = 0
	histo = 0.0
	ncontacts = 0
	
	! XXXX
	! XXXX Main routine....
	! XXXX

100	nframes=0
	framesdone = 0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	if (nframes.le.framestodiscard) goto 101

	framesdone = framesdone + 1

	! Loop over molecules of sp1
	i = s_start(sp1)-1
	do m1=1,s_nmols(sp1)

	  ! Loop over defined XH
	  do n=1,nxh

	    ! Get position of H atom
	    rh(1) = xpos(i+xh(n,2))
	    rh(2) = ypos(i+xh(n,2))
	    rh(3) = zpos(i+xh(n,2))

	    ! Get minimum image position of X atom, and H->X vector
	    call pbc(xpos(i+xh(n,1)),ypos(i+xh(n,1)),zpos(i+xh(n,1)),rh(1),rh(2),rh(3),rx(1),rx(2),rx(3))
	    vhx = rx - rh
	    dhx = sqrt(sum(vhx*vhx))
	!write(0,*) vhx,dhx
	    vhx = vhx / dhx

	    ! Loop over molecules of sp2
	    j = s_start(sp2)-1
	    do m2=1,s_nmols(sp2)

	      ! Loop over defined Y
	      do m=1,ny

		! Get minimum image position of Y w.r.t. H
		call pbc(xpos(j+y(m)),ypos(j+y(m)),zpos(j+y(m)),rh(1),rh(2),rh(3),ry(1),ry(2),ry(3))
		vhy = ry - rh
		dhy = sqrt(sum(vhy*vhy))
	!write(0,*) "Y",vhy,dhy
		vhy = vhy / dhy

		! Calculate XY distance
		dxy = dsqrt(sum((rx-ry)*(rx-ry)))
		
		! Calculate angle
		dp = sum(vhx*vhy)
		angle = safeAngle(dp)

		! Geometry check - is it a hydrogen bond by our definition?
		if ((dxy.gt.maxdist).or.(angle.lt.minang)) then
		  ! No longer an H-bond so, if flag was already set, add to histo
		  if (hbflag(m1,n,m2,m).ne.0) then
		    ncontacts = ncontacts + 1
		    histo(hbflag(m1,n,m2,m)) = histo(hbflag(m1,n,m2,m)) + 1.0
		    if (hbflag(m1,n,m2,m).gt.maxlife) maxlife = hbflag(m1,n,m2,m)
		    hbflag(m1,n,m2,m) = 0
		  end if
		  cycle
		end if

		! Increase flag counter for this contact
		hbflag(m1,n,m2,m) = hbflag(m1,n,m2,m) + 1

	      end do   ! Loop over Y

	      j = j + s_natoms(sp2)

	    end do   ! Loop over sp2 molecules

	  end do   ! Loop over XH

	  i = i + s_natoms(sp1)
	
	end do   ! Loop over sp1 molecules

	if (framesdone.EQ.framestodo) goto 800
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

	! Accumulate remaining H-bonds...
	i = 0
	do m1=1,s_nmols(sp1)
	  do n=1,nxh
	    do m2=1,s_nmols(sp2)
	      do m=1,ny
		if (hbflag(m1,n,m2,m).ne.0) then
		  i = i + 1
		  ncontacts = ncontacts + 1
		  histo(hbflag(m1,n,m2,m)) = histo(hbflag(m1,n,m2,m)) + 1.0
		  if (hbflag(m1,n,m2,m).gt.maxlife) maxlife = hbflag(m1,n,m2,m)
		  hbflag(m1,n,m2,m) = 0
		end if
	      end do   ! Loop over Y
	    end do   ! Loop over sp2 molecules
	  end do   ! Loop over XH
	end do   ! Loop over sp1 molecules
	if (i.gt.0) write(6,"(a,i8,a)") "Added on lifetime data for ", i, " contacts present at end of trajectory."

	! Write out information
	do n=1,maxlife
	  write(6,*) n, histo(n), histo(n) / ncontacts
	end do

999	write(0,*) "Finished."
	end program lifehb

