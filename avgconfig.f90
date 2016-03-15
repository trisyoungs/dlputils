!	** avgconfig.f90 **
!	Calculate the average coordinates over a series of trajectory frames

	program avgconfig
	use dlprw; use utility
	implicit none
	integer :: maxframes,vec(3),baselen,nargs,success,n,i,targetsp
	integer :: nframes, nframesused, xbin, ybin, zspecies
	character*80 :: hisfile,outfile,basename,resfile,temp,flagfile,altheaderfile
	character*8 :: discard
	logical :: altheader = .FALSE.
	real*8, allocatable :: pos(:,:,:), avg(:,:), sd(:,:), minimum(:,:), maximum(:,:)
	real*8 :: tx, ty, tz, total(3)
	integer :: iargc

	write(0,*) "*** avgconfig"

	nargs = iargc()
	if (nargs.LT.2) then
	  write(*,"(a)") "Usage: avgconfig <HISTORYfile> <OUTPUTfile> [-options]"
	  write(*,"(a)") "        [-frames n]          Number of frames to do"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
	end if

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  IF (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 50
	  endIF
	end do
50	if (baselen.eq.-1) then
	  basename="pointresults."
	  baselen=14
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=20,file=basename(1:baselen)//".out",status='replace',form='formatted')

	write(0,"(A,A)") "History file : ",hisfile
	write(20,"(A,A)") "History file : ",hisfile
	!call openhis(hisfile,n)
	write(0,"(A,A)") " Output file : ",outfile
	write(20,"(A,A)") " Output file : ",outfile
	if (outinfo(outfile,1).eq.-1) goto 798

	! Set some variable defaults before we read in any command-line arguments
	maxframes = 1000

	n = 2
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-frames") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") maxframes
	      write(0,"(A,I1)") "Frames to do is", maxframes
	    case default
	      write(0,*) "Unrecognised command line option : ",temp
	      stop
	  end select
	end do
	
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

	! Allocate arrays
	allocate(pos(maxframes,natms,3))
	allocate(avg(natms,3))
	allocate(sd(natms,3))
	allocate(minimum(natms,3))
	allocate(maximum(natms,3))

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes = 0

102	success=readframe()
	IF (success.eq.1) goto 801  ! End of file encountered....
	IF (success.eq.-1) goto 119  ! File error....
	nframes=nframes+1
	if (MOD(nframes,100).eq.0) write(0,*) nframes
	
	! If this is the first frame, store the position.
	! Otherwise, store MIM position
	if (nframes.eq.1) then
	  ! Store positions from this frame
	  pos(nframes,:,1) = xpos(:)
	  pos(nframes,:,2) = ypos(:)
	  pos(nframes,:,3) = zpos(:)
	else
	  do i=1,natms
            call pbc(xpos(i),ypos(i),zpos(i), pos(1,i,1), pos(1,i,2), pos(1,i,3), tx, ty, tz)
	    pos(nframes,i,1) = tx
	    pos(nframes,i,2) = ty
	    pos(nframes,i,3) = tz
	  end do
	end if

	! Next frame (or finish)
116	if (nframes.EQ.maxframes) goto 801 
	goto 102

119	write(0,*) "HISTORY file ended prematurely!"
	write(0,"(A,I5,A)") "Managed ",nframes," frames before error."

	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "Problem with HISTORY file - failed to read header."
	goto 999
800	write(0,*) "End of unformatted HISTORY file found."
801	write(0,*) "Finished."

	! Take averages and standard deviation
	do i=1,natms
          minimum(i,1) = minval(pos(:,i,1))
          maximum(i,1) = maxval(pos(:,i,1))
          minimum(i,2) = minval(pos(:,i,2))
          maximum(i,2) = maxval(pos(:,i,2))
          minimum(i,3) = minval(pos(:,i,3))
          maximum(i,3) = maxval(pos(:,i,3))
          avg(i,1) = sum(pos(:,i,1)) / real(nframes)
          avg(i,2) = sum(pos(:,i,2)) / real(nframes)
          avg(i,3) = sum(pos(:,i,3)) / real(nframes)
          ! Now for S.D.
          total = 0.0d0
          do n=1,nframes
            total(1) = total(1) + (pos(n,i,1) - avg(i,1))**2
            total(2) = total(2) + (pos(n,i,2) - avg(i,2))**2
            total(3) = total(3) + (pos(n,i,3) - avg(i,3))**2
          end do
          sd(i,1) = sqrt( total(1) / real(nframes) )
          sd(i,2) = sqrt( total(2) / real(nframes) )
          sd(i,3) = sqrt( total(3) / real(nframes) )
	end do

	! Write out data
	open(unit=11,form="formatted",file=basename(1:baselen)//"avg.CONFIG",status="replace")
	write(11,"(a,i6,a)") "Average of ",nframes," configurations"
	write(11,"(2i10)") 0, imcon
	if (imcon.gt.0) then
	  write(11,"(3f20.14)") cell
	end if
	do i=1,natms
	  write(11,"(a8)") atmname(i)
	  write(11,"(3f20.14)") avg(i,:)
	end do

	open(unit=12,form="formatted",file=basename(1:baselen)//"avg.akf",status="replace")
	if (imcon.gt.0) then
	  write(12,"('cellmat  ',9f12.6)") cell
	end if
	do i=1,natms
	  write(12,"('atom  ',a8,2x,3f12.6)") atmname(i), avg(i,:)
	  write(12,"('sphere  ',i6,3f12.6)") i, sd(i,:)
	end do

	close(11)
	close(12)

	write(0,*) "Finished!"
999	close(10)
	close(20)

	end program avgconfig
