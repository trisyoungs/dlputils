!	** cn **
!	Calculate CNs based on distance between specified species com's or atom centres

	program cn
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,a1,s1,sp1,sp2,m1,m2,baselen,nframes,success,nargs,divide = -1, compairs(10,2)
	integer :: count, start, finish, framestodo = -1, framesdone, frameskip = 0
	integer :: iargc
	logical :: altheader = .FALSE.
	real*8 :: c1x,c1y,c1z,c2x,c2y,c2z,tx,ty,tz,rij,cutoff, total,raverage

	nargs = iargc()
	if (nargs.lt.5) then
	  write(0,"(A)") "Usage : cn <DLP HISTORYfile> <DLP OUTPUTfile> <sp1> <sp2> <cutoff> [-header file] [-frames n] [-skip n] [-divide n] [-compair sp i j]"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(I4)") sp1
	call getarg(4,temp); read(temp,"(I4)") sp2
	call getarg(5,temp); read(temp,"(F10.8)") cutoff

	n = 5
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
	      altheader = .TRUE.
            case ("-frames")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
              write(0,"(A,I4)") "Frames to do: ",framestodo
            case ("-skip")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") frameskip
              write(0,"(A,I4)") "Frames to skip at start: ",frameskip
            case ("-divide")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") divide
              write(0,"(A,I4)") "Output will be partitioned into chunk of (frames) ",divide
            case ("-compair")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") s1
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(s1,1)
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(s1,2)
              write(0,"(A,3I4)") "Using COMpair for species ",s1, compairs(s1,:)
	    case default
	      write(0,"(a,a)") "Unrecognised command line option:",temp
	      stop
	  end select
	end do

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) then
	  if (altheader) then
            write(0,*) "Restarted trajectory:"
            close(dlpun_his)
            call openhis(altheaderfile,10)
            if (readheader().EQ.-1) goto 799
            close(dlpun_his)
            call openhis(hisfile,10)
          else
	    goto 799
	  end if
	end if

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

	resfile=basename(1:baselen)//"cn"//CHAR(48+sp1)//CHAR(48+sp2)
	open(unit=9,file=resfile,form="formatted",status="replace")
	write(9,"(a,i2,a,i2,a,f10.4)") "# CN of species ",sp2," around ",sp1," within ",cutoff

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
	raverage = 0.0
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
	do s1=1,nspecies
	  if (compairs(s1,1).ne.0) then
	    do m1=1,s_nmols(s1)
	      c1x = xpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,1)-1)
	      c1y = ypos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,1)-1)
	      c1z = zpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,1)-1)
	      c2x = xpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,2)-1)
	      c2y = ypos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,2)-1)
	      c2z = zpos(s_start(s1)+(m1-1)*s_natoms(s1)+compairs(s1,2)-1)
	      call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	      comx(s1,m1) = (c1x+tx)*0.5
	      comy(s1,m1) = (c1y+ty)*0.5
	      comz(s1,m1) = (c1z+tz)*0.5
	    end do
	  end if
	end do

	total = 0.0
	do m1=1,s_nmols(sp1)
	  ! Get the centre of mass for species sp1 mol m1...
	  c1x=comx(sp1,m1)
	  c1y=comy(sp1,m1)
	  c1z=comz(sp1,m1)
	  do m2=1,s_nmols(sp2)
	    ! Must do a quick check to exclude m1=m2 if sp1=sp2
	    if ((sp1.EQ.sp2).AND.(m1.EQ.m2)) cycle
	    c2x=comx(sp2,m2)
	    c2y=comy(sp2,m2)
	    c2z=comz(sp2,m2)
	    ! Get the shortest (MIM) distance between the pair of positions...
	    call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	    rij=sqrt( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )
	    ! Store this distance...
	    if (rij.lt.cutoff) total = total + 1.0
	    !if (rij.lt.cutoff) write(0,*) nframes,m2, rij
	    !if (rij.lt.cutoff) write(0,"(6f12.6)") c1x,c1y,c1z,tx,ty,tz
	  end do
	end do

	! Normalise
	if (sp1.EQ.sp2) then
	  total = total / (s_nmols(sp1)-1)
	else
	  total = total / s_nmols(sp1)
	end if
	
	! Increase total accumulator
	raverage = raverage + total

	if (divide.eq.-1) then
	  write(9,"(2F10.6)") total,raverage/framesdone
	else if (mod(framesdone,divide).eq.0) then
	  write(9,"(i8,f10.6)") framesdone, raverage/divide
	  raverage = 0.0
	end if

	if (framesdone.EQ.framestodo) goto 800
	! Next frame
	goto 101

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended before framestodo was fulfilled..."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",framesdone," (wanted ",framestodo,")"
	goto 801
800	write(0,*) "Framestodo was fulfilled."
801	write(0,*) ""

	close(9)
	write(0,*) "Finished."
999	close(10)
	close(13)
	end program cn
