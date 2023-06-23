	! #####################
	! XYZ to HISu converter
	! #####################

	program xyz2his
	use dlprw; use parse; use ptable
	implicit none
	character*80 :: xyzfile,temp,outfile,header
	real*8 :: a
	integer :: nargs,n,i,m,success,nframes,destfmt
	integer :: iargc, nAssignedMasses = 0
	logical :: hascell = .FALSE., extraline = .FALSE., found
	character*8, allocatable :: assignedMasses(:)

	nargs = iargc()
	if (nargs.lt.2) stop " Usage: xyz2his <xyz file> <output.HISu> [-cell ax ay az bx by bz cx cy cz] [-extraline] [-cubic a]"
	call getarg(1,xyzfile)
	call getarg(2,outfile)
	imcon = 0

	n = 2
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-cell")
	      do m=1,9
                n = n + 1
	        if (n.GT.nargs) stop "Error: Partial cell specification given?"
	        call getarg(n,temp); read(temp, "(f20.14)") cell(m)
	      end do
	      hascell = .TRUE.
	      imcon = 3
              write(0,"(A)") "Read full cell specification as:"
	      write(0,"(3f12.4)") cell
	    case ("-cubic")
	      n = n + 1; 
	      call getarg(n,temp); read(temp, "(f20.14)") a
	      cell = 0.0
	      cell(1) = a
	      cell(5) = a
	      cell(9) = a
	      hascell = .TRUE.
	      imcon = 1
	      write(0,*) "Cubic cell will be written, with side length ", cell(1)
            case ("-extraline")
	      extraline = .TRUE.
              write(0,"(A)") "Assuming presence of blank line in between frames"
	    case default
	      write(0,"(a,a)") "Unrecognised command line option:",temp
	      stop
	  end select
	end do

	! Set destination format (0 = unformatted history file, 1 = formatted history file)
	destfmt = 0
	
	! Open the xyz file and read in the first frame to get atom numbers and names
	open(unit=15,file=xyzfile,form='formatted',status='old',err=999)
	read(15,"(i20)") natms
	call alloc_xyz
	read(15,"(a80)") header
	do n=1,natms
	  if (.not.readline(15)) goto 999
	  atmname(n) = arg(1)
	  keytrj = 0
	  if (nargsparsed.ge.7) keytrj = 1
	  if (nargsparsed.eq.10) keytrj = 2
	end do
	rewind(15)

	! Handle atomic masses
	allocate(assignedMasses(1:natms))
	do n=1,natms
	  ! Assign mass for the current atom name
	  mass(n) = getMass(atmname(n))
	  ! If we haven't already assigned this atom name, add it to our array
	  found = .FALSE.
	  do m=1,nAssignedMasses
	    if (assignedMasses(m).eq.atmname(n)) then
	      found = .TRUE.
	      exit
	    end if
	  end do
	  if (.not.found) then
	    nAssignedMasses = nAssignedMasses + 1
	    assignedMasses(nAssignedMasses) = atmname(n)
	    write(6,*) " -- Atom name ", atmname(n), "assigned mass of ", mass(n)
	  end if
	end do

	nframes = 0
	write(0,*) "Writing header of new file..."
	if (writeheader(outfile,14,destfmt).EQ.-1) stop "Failed to write new history file header!"

10	read(15,"(i20)",end=15,err=15) natms

	! Check number of atoms - if zero, we probably tried to read something other than the natms integer
	if (natms.eq.0) then
	  write(0,*) "Error reading natoms for frame ", nframes+1
	  goto 15
	end if
	read(15,"(a80)") header
	do n=1,natms
	  if (.not.readline(15)) then
	    write(0,*) "Error reading frame ", nframes+1, " at atom ", n, "(premature end of file?)"
	    goto 999
	  end if
	  if (nargsparsed.lt.4) then
	    write(0,*) "Error reading frame ", nframes+1, " at atom ", n, "(less than three numbers follow the atom name)"
	    write(0,*) "Number of args read = ", nargsparsed
	    write(0,"(i4,2x,a80)") (m, arg(m), m=1,nargsparsed)
	    goto 999
	  end if
	  xpos(n) = argr(2)
	  ypos(n) = argr(3)
	  zpos(n) = argr(4)
	  if (nargsparsed.ge.7) then
	    xvel(n) = argr(5)
	    yvel(n) = argr(6)
	    zvel(n) = argr(7)
	  end if
	  if (nargsparsed.eq.10) then
	    xfor(n) = argr(8)
	    yfor(n) = argr(9)
	    zfor(n) = argr(10)
	  end if
	end do

	nframes = nframes + 1
	if (MOD(nframes,100).EQ.0) write(0,*) nframes

	! Write out the frame...
	success = writeframe()
	if (success.EQ.-1) stop "Failed to write frame!"

	! Skip extra blank line if requested
	if (extraline) then
	  if (.NOT.readline(15)) goto 15
	end if

	goto 10
15	close(15)  ! end of file..
	write(0,"(A,I5,A,A)") "Read/wrote ",nframes," from file ",xyzfile
	  
20      close(14)
	stop "Finished conversion."
999	stop "Something bad happened!"

	end program xyz2his
