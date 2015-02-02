!	** surfacify **
!	Calculate surface 'density' of atoms per frame of simulation box

	program surfacify
	use dlprw; use utility
	implicit none
	integer, parameter :: MAXSPECIES = 10
	real*8, allocatable :: pdens(:,:,:)
	integer :: atoms(MAXSPECIES,100), grid(3)
	real*8 :: delta(3), tx, ty, tz
	integer :: nmapframes, baselen,nargs,nframesused, collect = 1, framescollected
	character*80 :: hisfile,outfile,basename,resfile,temp,flagfile,altheaderfile
	logical :: altheader = .FALSE.
	integer :: success,n1,n2,n3, i, n, nframes, aoff1
	integer :: sp1,m1,startf,endf
	integer :: iargc,tth,th,hun,ten,units

	write(0,*) "*** surfacify"

	nargs = iargc()
	if (nargs.LT.2) then
	  write(*,"(a)") "Usage: surfacify <HISTORYfile> <OUTPUTfile> [options...]"
	  write(*,"(a)") "        [-atom sp i]            Include atom i of species sp in the surface"
	  write(*,"(a)") "        [-grid x y z]           Grid to use in each direction"
	  write(*,"(a)") "        [-start frameno]        Trajectory frame to start calculations (default = 1)"
	  write(*,"(a)") "        [-end frameno]          Trajectory frame to end calculations on (default = last)"
	  write(*,"(a)") "        [-collect nframes]      Sum distribution over clusters of consecutive frames"
	  write(*,"(a)") "        [-header file]          Use specified file to get header"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
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
	  basename="3ddistresults."
	  baselen=14
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=15,file=basename(1:baselen)//"surfacify",form="formatted",status="replace")

	startf = 1
	endf = 0
	atoms = 0
	collect = 1

	n = 2
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-atom") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp1
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") i
	      atoms(sp1,i) = 1
	      write(0,"(A,I3,A,I2)") "Added atom",i," from species ",sp1
	    case ("-grid")
	      n = n + 1; call getarg(n,temp); read(temp,"(i6)") grid(1)
	      n = n + 1; call getarg(n,temp); read(temp,"(i6)") grid(2)
	      n = n + 1; call getarg(n,temp); read(temp,"(i6)") grid(3)
	      write(0,"(A,3i5)") "Grid points = ",grid
	    case ("-start")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") startf
	      write(0,"(A,I5)") "Starting frame = ",startf
	    case ("-end")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") endf
	      write(0,"(A,I5)") "End frame = ",endf
	    case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
              altheader = .TRUE.
	    case ("-collect")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") collect
	      write(0,"(A,I5)") "Density will be be taken over framesets of size = ",collect
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

	! Print out a summary of the control variables to the output file.
	write(15,"(A,A)") "Input file: ",hisfile
	write(15,"(A,I3)") "Molecular species in file : ",nspecies
	do sp1=1,nspecies
	  ! Check that molecules have had their axes defined
	  write(15,"(A,i3)") "...atoms for species ",sp1
	  write(15,"(50i3)") atoms(sp1,:)
	end do
	delta(1) = cell(1) / grid(1)
	delta(2) = cell(5) / grid(2)
	delta(3) = cell(9) / grid(3)
	write(15,"(A,3F6.3)") "Grid spacings = ",delta
	write(15,"(A,3I4)") "Grid (intermolecular) points in each XYZ = ",grid
	write(15,"(A,I5,A,I5,A)") "Frame range = ",startf," to ",endf," (0=last)"

	allocate(pdens(0:grid(1),0:grid(2),0:grid(3)))

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
	nframes = 0
	nframesused = 0
	pdens = 0.0
	framescollected = 0

101	success=readframe()
	IF (success.eq.1) goto 120  ! End of file encountered....
	IF (success.lt.0) goto 119  ! File error....
	nframes=nframes+1

	if (MOD(nframes,100).eq.0) then
	  write(0,"(i6)") nframes
	end if
	if (nframes.LT.startf) goto 101

	nframesused = nframesused + 1

	aoff1 = 0
	do sp1=1,nspecies
	  do m1=1,s_nmols(sp1)
	    do i=1,s_natoms(sp1)

	      if (atoms(sp1,i).eq.0) cycle

	      tx = xpos(aoff1+i)
	      ty = ypos(aoff1+i)
	      tz = zpos(aoff1+i)

	      call unitfold(tx,ty,tz)

	      ! Calculate integer position in grid
	      n1=(tx/delta(1))+1
	      n2=(ty/delta(2))+1
	      n3=(tz/delta(3))+1

	      pdens(n1,n2,n3) = pdens(n1,n2,n3) + 1

	    end do
	    aoff1 = aoff1 + s_natoms(sp1)
	  end do
	end do

	framescollected = framescollected + 1

	! Ready to save density?
	if (framescollected.eq.collect) then

	  ! Normalise
	  pdens(n1,n2,n3)=pdens(n1,n2,n3)/(delta(1)*delta(2)*delta(3)) / (collect*1.0)

	  tth = nframes / 10000; n = nframes - tth*10000
	  th = n / 1000; n = n - th*1000
	  hun = n / 100; n = n - hun*100
	  ten = n / 10; n = n - ten*10
	  units = n
	  resfile=basename(1:baselen)//CHAR(48+tth)//CHAR(48+th)//CHAR(48+hun)//CHAR(48+ten)//CHAR(48+units)//"_sfcfy.pdens"
	  open(unit=9,file=resfile,form="formatted")
	  write(9,*) grid(1)+1,grid(2)+1,grid(3)+1
	  write(9,"(9f8.4)") delta(1),0.0,0.0,0.0,delta(2),0.0,0.0,0.0,delta(3)
! 	  write(9,"(3f10.4)") -cell(1)/2.0,-cell(5)/2.0,-cell(9)/2.0
 	  write(9,"(3f10.4)") 0.0,0.0,0.0
	  write(9,*) "zyx"
	  do n1=0,grid(1)
	    do n2=0,grid(2)
	      do n3=0,grid(3)
	  	write(9,"(f12.8)") pdens(n1,n2,n3)
	      end do
	    end do
	  end do
	  close(9)

	  framescollected = 0
	  pdens = 0.0
	end if

	! Next frame (or finish)
116	if (nframes.EQ.endf) goto 120
	goto 101

119	write(0,*) "HISTORY file ended prematurely!"
	write(0,"(A,I5,A)") "Managed ",nframesused," frames before error."
	write(15,"(A)") "HISTORY file ended prematurely!"
	write(15,"(A,I5,A)") "Managed ",nframesused," frames before error."
120	write(0,*) "Finished."
	write(15,"(A)") "Finished."
	write(15,"(A,I5,A)") "Averages will be taken over ",nframesused," frames."

	goto 999

	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "Problem with HISTORY file - failed to read header."
	goto 999
	write(0,*) "End of unformatted HISTORY file found."
	write(15,"(A)") "End of unformatted HISTORY file found."

	write(0,*) "Finished!"
	write(15,"(A)") "Finished!"
999	close(10)
	close(13)

	end program surfacify
