	! Parser Module

	module parse

	  integer, parameter :: MAXARGS = 50, MAXARGLEN = 30
	  integer :: nargsparsed
	  character(MAXARGLEN) :: arg(MAXARGS)

	contains

	! Read line from specified unit number and split into delimited arguments
	logical function readline(uno)
	implicit none
	character*200 line
	integer n,uno,arglen
	logical done
	! Blank previous arguments
	do n=1,MAXARGS
	  arg(1:MAXARGLEN) = " "
	end do
	! Read line
	read(uno,"(a200)",end=900,err=900) line
	! Go through characters in line and get arguments
	nargsparsed = 0
	arglen = 0
	done = .false.
	do n=1,200
	  select case (ichar(line(n:n)))
	    case (32,44,9)
	      ! Blank character (space=32, comma=44, tab=9)
	      ! If arglen != 0 then 'store' current argument.
	      ! Otherwise, move on
	      if (arglen.ne.0) then
		arglen = arglen + 1
		! arg(nargsparsed+1)(arglen:arglen) = ENDOFLINE
		nargsparsed = nargsparsed + 1
		arglen = 0
	      end if
	    case (33,35)
	      ! Comment character - (bang=33, hash=35) we're done
	      done = .true.
	    case default
	      ! Add character to argument
	      arglen = arglen + 1
	      arg(nargsparsed+1)(arglen:arglen) = line(n:n)
	  end select
	  if (done) exit
	end do
	goto 999
900	nargsparsed = -1
	readline = .FALSE.
	return
999	readline = .TRUE.
	end function readline

	! Return argument as integer
	integer function argi(n)
	implicit none
	integer n
	if (n.gt.MAXARGS) then
	  write(0,*) "argi() was passed an argument number greated than MAXARGS - ",n
	  argi = 0
	else
	  read(arg(n),"(i20)") argi
	end if
	end function argi

	! Return argument as real*8
	real*8 function argr(n)
	implicit none
	integer n
	if (n.gt.MAXARGS) then
	  write(0,*) "argd() was passed an argument number greated than MAXARGS - ",n
	  argr = 0.0
	else
	  read(arg(n),"(f30.10)") argr
	end if
	end function argr

	! Retrieve CLI argument specified as a floating-point number
	real*8 function getargr(argindex)
	implicit none
	integer, intent(in) :: argindex
	character*80 :: temp
	integer :: nargs, iconvert

	getargr = 0.0

	! Check argument index, and get argument if OK
	if (argindex.gt.iargc()) then
	  write(0,*) "Argument index", argindex, "out of range in getargr."
	  stop
	endif
	call getarg(argindex,temp)

	! Check contents, and convert
	if ((index(temp,'.').gt.0).or.(index(temp,'e').gt.0).or.(index(temp,'E').gt.0)) then
	  read(temp,*) getargr
	else
	  read(temp,*) iconvert
	  getargr = iconvert
	end if

	end function getargr

	! Retrieve CLI argument specified as an integer number
	real*8 function getargi(argindex)
	implicit none
	integer, intent(in) :: argindex
	character*80 :: temp
	integer :: nargs
	real*8 :: rconvert

	getargi = 0.0

	! Check argument index, and get argument if OK
	if (argindex.gt.iargc()) then
	  write(0,*) "Argument index", argindex, "out of range in getargi."
	  stop
	endif
	call getarg(argindex,temp)

	! Check contents, and convert
	if ((index(temp,'.').gt.0).or.(index(temp,'e').gt.0).or.(index(temp,'E').gt.0)) then
	  read(temp,*) rconvert
	  getargi = int(rconvert)
	else
	  read(temp,*) getargi
	end if

	end function getargi

	end module parse
