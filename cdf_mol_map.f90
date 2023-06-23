!	** cdf_mol_map **
!       Generates a molmap file (using -1 to indicate end of frame) by calculation
!	cylindrical distribution functions between the centres-of-mass of
!	molecules and a specified vector and using max and min distance to only select molecules in that "ring"

	program cyldf_molmap
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp, mapfile
	logical :: altheader = .FALSE.
	integer :: nbins,n,m,a1,s1,m1,baselen,bin,nframes,success,nargs,numadded,framestodo = -1,framestodiscard = 0,framesdone, compairs(10,2),sp
	integer :: iargc
	real*8 :: dist,pos(3),integral,origin(3),vector(3),t1(3),t2(3)
	real*8 :: maxdist, mindist, xyyx, xzzx, yzzy, denom, numdens, shellvol
	real*8, allocatable :: cdf(:,:)

	!binwidth=0.1   ! In Angstroms
	!compairs = 0
	nargs = iargc()
	if (nargs.LT.8) stop "Usage : cdf_mol_map <HISTORYfile> <OUTPUTfile> <ox> <oy> <oz> <vx> <vy> <vz> <sp> [-maxdist maxdist] [-mindist mindist] [-header hisfile] [-frames n] [-discard n] [-mapfile mapfile]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(f20.14)") origin(1)
	call getarg(4,temp); read(temp,"(f20.14)") origin(2)
	call getarg(5,temp); read(temp,"(f20.14)") origin(3)
	call getarg(6,temp); read(temp,"(f20.14)") vector(1)
	call getarg(7,temp); read(temp,"(f20.14)") vector(2)
	call getarg(8,temp); read(temp,"(f20.14)") vector(3)
        call getarg(9,temp); read(temp,"(i4)") sp
	write(0,"(a,3f10.4)") "Vector origin is : ", origin
	write(0,"(a,3f10.4)") "Vector is : ", vector
	write(0,"(a,i4)") "Species is :", sp

	! Check line parameters
	denom = dot_product(vector,vector)
	if (denom.lt.1.0e-6) stop "Invalid line vector supplied."

        !Set default values
        mapfile="map.out"
        mindist=0.0
        maxdist=1000.0
	
	n = 9
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-maxdist")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") maxdist
	      write(0,"(A,f6.2)") "Max distance for molmap set to ",maxdist
            case ("-mindist")
	      n = n + 1; call getarg(n,temp); read(temp,"(F10.4)") mindist
	      write(0,"(A,f6.2)") "Min distance for molmap set to ",mindist
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
            case ("-mapfile")
              n = n + 1; call getarg(n,temp); read(temp,"(A)") mapfile
            case default
	      write(0,"(a,a)") "Unrecognised command line option:",temp
	      stop
	  end select
	end do

        write(0,"(A)") "Output mol_map file set to ",mapfile

        !check that a min or max distance has been specified
        if((mindist.lt.0.01).AND.(maxdist.gt.999.0)) then
        stop "no max or min distance specified"
        endif
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

        !open map file for writting
	open(unit=9,file=mapfile,form="formatted")
 
	!nbins = maxval(cell) / binwidth + 1
	!write(0,"(A,I5,A,F6.3,A)") "There will be ",nbins," histogram bins of ",binwidth," Angstroms."
	!allocate(cdf(nspecies,nbins))
	!cdf = 0.0

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
100	nframes=0
	framesdone = 0
101     success=readframe()
         !write(0,*) "read frame"
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes, framesdone
	if (nframes.le.framestodiscard) goto 101

	framesdone = framesdone + 1

	call calc_com
        !write(0,*) "got to here"

	  write(9,"(i5,a)",advance='no') framesdone, " "
	  do m1=1,s_nmols(sp)     ! Loop over all molecules of species sp
	    ! call pbc(comx(s1,m1),comx(s1,m1),comy(s1,m1),origin(1),origin(2),origin(3),t1(1),t1(2),t1(3))
	    pos(1) = comx(sp,m1) - origin(1)
	    pos(2) = comy(sp,m1) - origin(2)
	    pos(3) = comz(sp,m1) - origin(3)
	    xyyx = vector(1)*pos(2) - vector(2)*pos(1)
	    xzzx = vector(1)*pos(3) - vector(3)*pos(1)
	    yzzy = vector(2)*pos(3) - vector(3)*pos(2)
	    t1(1) = vector(2)*xyyx + vector(3)*xzzx
	    t1(2) = vector(3)*yzzy - vector(1)*xyyx
	    t1(3) = -vector(1)*xzzx - vector(2)*yzzy
	    dist = sqrt(sum(t1*t1)) / denom
            
            !write(0,*) "
            if ((dist.GT.mindist).AND.(dist.LT.maxdist)) then
 	    write (9,"(i5,a)",advance='no') m1, " "
	    endif

	  end do
        write (9,"(i5)") -1
	!end do

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

	dist = vector(1)*cell(1) + vector(2)*cell(5) + vector(3)*cell(9)
	
        write(0,*) "Finished."
999	close(10)
	close(13)
        close(9)
	end program cyldf_molmap

