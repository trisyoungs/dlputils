!	** zangle.f90 **
!	Calculate angle between XY planes (from axis definition)

	program zangle
	use dlprw; use utility
	implicit none
	real*8 :: ahist(0:180)
	integer :: baselen,nargs,nframesused
	character*80 :: hisfile,outfile,basename,resfile,temp,altheaderfile
	character*8 :: discard
	logical :: altheader = .FALSE.
	integer :: success, bin, n, m, nframes
	real*8 :: numadded, angle, mindist, maxdist, dist, tx,ty,tz, sumy, sumy2
	integer :: sp1,m1,sp2,m2,startf,endf
	integer :: iargc

	write(0,*) "*** zangle"

	nargs = iargc()
	if (nargs.lt.4) then
	  write(*,"(a)") "Usage: zangle <HISTORYfile> <OUTPUTfile> <sp1> <sp2>"
	  write(*,"(a)") "        [-axis sp x1 x2 y1 y2]  Atoms to use for axis calculation in species sp"
	  write(*,"(a)") "        [-start frameno]        Trajectory frame to start calculations (default = 1)"
	  write(*,"(a)") "        [-end frameno]          Trajectory frame to end calculations on (default = last)"
	  write(*,"(a)") "        [-mindist r]            Minimum separation allowed between molecules (default = 0.0)"
	  write(*,"(a)") "        [-maxdist r]            Maximum separation allowed between molecules (default = 10000.0)"
	  write(*,"(a)") "        [-header file]          Use specified file to get header"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
	  call getarg(3,temp); read(temp,"(I4)") sp1
	  call getarg(4,temp); read(temp,"(I4)") sp2
	end if
	n = 10
	write(0,"(A,A)") "History file : ",hisfile
	!call openhis(hisfile,n)
	write(0,"(A,A)") " Output file : ",outfile
	if (outinfo(outfile,1).eq.-1) goto 798

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  IF (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 50
	  endIF
	end do
50	if (baselen.eq.-1) then
	  basename="adensresults."
	  baselen=13
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=15,file=basename(1:baselen)//"out",form="formatted",status="replace")

	! Set some variable defaults before we read in any command-line arguments
	call alloc_axis
	mindist = 0.0
	maxdist = 10000.0
	startf = 1
	endf = 0
	numadded = 0

	n = 4
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-axis") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp1
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp1,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp1,2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp1,3)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") axesAatoms(sp1,4)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",sp1," calculated from: X=",axesAatoms(sp1,1),"->", &
	        & axesAatoms(sp1,2),", Y=0.5X->0.5(r(",axesAatoms(sp1,3),")->r(",axesAatoms(sp1,4),"))"
	      do m=1,3
		if ((axesAatoms(sp1,m).lt.1).or.(axesAatoms(sp1,m).gt.s_natoms(sp1))) stop "Atom id out of range for axes on this species!"
	      end do
	      axesAdefined(sp1) = .true.
	    case ("-start")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") startf
	      write(0,"(A,I5)") "Starting frame = ",startf
	    case ("-end")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") endf
	      write(0,"(A,I5)") "End frame = ",endf
	    case ("-mindist")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") mindist
	      write(0,"(A,f8.4)") "Minimum separation between molecules = ",mindist
	    case ("-maxdist")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") maxdist
	      write(0,"(A,f8.4)") "Maximum separation between molecules = ",maxdist
	    case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
              altheader = .TRUE.
	    case default
	      write(0,*) "Unrecognised command line option : ",temp
	      stop
	  end select
	end do
	
	! Print out a summary of the control variables to the output file.
	write(15,"(A,A)") "Input file: ",hisfile
	write(15,"(A,I3)") "Molecular species in file : ",nspecies
	do n=1,nspecies
	  ! Check that molecules have had their axes defined
	  if (axesAdefined(n)) then
	    write(15,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Local axes for species ",n," calculated from: X=",axesAatoms(n,1),"->", &
	        & axesAatoms(n,2),", Y=0.5X->0.5(r(",axesAatoms(n,3),")->r(",axesAatoms(n,4),"))"
	  else
	    if ((n.eq.sp1).or.(n.eq.sp2)) stop "Axes must be defined for the two target species."
	  end if
	end do
	write(15,"(A,I4)") "Central species = ",sp1
	write(15,"(A,I4)") "Surrounding species = ",sp2
	write(15,"(A,I5,A,I5,A)") "Frame range = ",startf," to ",endf," (0=last)"

	! Clear the arrays before we start...
	write(0,*) "Clearing initial arrays..."
	ahist = 0.0
	numadded = 0.0

	! Read the header of the history file...
        call openhis(hisfile,10)
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

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes = 0
	nframesused = 0
102	success=readframe()
	IF (success.eq.1) goto 120  ! End of file encountered....
	IF (success.eq.-1) goto 119  ! File error....
	nframes=nframes+1
	if (MOD(nframes,100).eq.0) write(0,"(i6)") nframes
	if (nframes.LT.startf) goto 102

	nframesused = nframesused + 1

	! Generate all molecular axes
	call genaxes()

	! Calculate histogram
	do m1=1,s_nmols(sp1)      ! Loop over all molecules of species 1...

	  do m2=1,s_nmols(sp2)

	    ! Don't consider same molecule with itself
	    if ((m1.eq.m2).and.(sp1.eq.sp2)) cycle

	    ! Get minimum image position of second centre with first
	    call pbc(axesAorigin(sp2,m2,1),axesAorigin(sp2,m2,2),axesAorigin(sp2,m2,3),axesAorigin(sp1,m1,1),axesAorigin(sp1,m1,2),axesAorigin(sp1,m1,3),tx,ty,tz)
	    tx = tx - axesAorigin(sp1,m1,1)
	    ty = ty - axesAorigin(sp1,m1,2)
	    tz = tz - axesAorigin(sp1,m1,3)
	    dist = sqrt(tx*tx + ty*ty + tz*tz)
	    if (dist.lt.mindist) cycle
	    if (dist.gt.maxdist) cycle

	    ! This molecule is within the limits defined, so get dot product of z-axes
	    angle = safeAngle( axesA(sp2,m2,7)*axesA(sp1,m1,7) + axesA(sp2,m2,8)*axesA(sp1,m1,8) + axesA(sp2,m2,9)*axesA(sp1,m1,9) )
	    bin = int(angle) 
	    ahist(bin) = ahist(bin) + 1.0

	    numadded = numadded + 1

	  end do

	end do

	! Next frame (or finish)
116	if (nframes.EQ.endf) goto 120
	goto 102

119	write(0,*) "HISTORY file ended prematurely!"
	write(0,"(A,I5,A)") "Managed ",nframesused," frames before error."
	write(15,"(A)") "HISTORY file ended prematurely!"
	write(15,"(A,I5,A)") "Managed ",nframesused," frames before error."
120	write(0,*) "Finished."
	write(15,"(A)") "Finished."
	write(15,"(A,I5,A)") "Averages will be taken over ",nframesused," frames."

	! ### Normalise histogram

	ahist = ahist / numadded

	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "Problem with HISTORY file - failed to read header."
	goto 999
800	write(0,*) "End of unformatted HISTORY file found."
	write(15,"(A)") "End of unformatted HISTORY file found."
801	write(0,"(A)") "-- Writing output files..."

	resfile=basename(1:baselen)//"zang"//CHAR(48+sp1)//CHAR(48+sp2)
	open(unit=11,file=resfile,form='formatted',status='replace')
	write(11,"('# Normalised histogram - caught ',i10,' molecules over ',i7,' frames')") numadded, nframesused
	write(11,"('# ZAngle calculated for species ',i2,' around species ',i2)") sp1, sp2
	write(11,"('# Distance min / max = ',f10.5,'/',f10.5)") mindist,maxdist
	sumy = 0.0
	sumy2 = 0.0
	do n=0,180
	  sumy = sumy + ahist(n)
	  sumy2 = sumy2 + ahist(n) + ahist(180-n)
	  if (n.eq.90) sumy2 = sumy2 - ahist(90)
	  if (n.le.90) then
	    write(11,"(5f10.4)") real(n), ahist(n), sumy, ahist(n) + ahist(180-n), sumy2
	  else
	    write(11,"(3f10.4)") real(n), ahist(n), sumy
	  end if
	end do

	write(0,*) "Finished!"
	write(15,"(A)") "Finished!"
999	close(10)
	close(11)

	end program zangle
