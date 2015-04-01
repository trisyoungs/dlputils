	! Parser Module

	module parse

	  integer, parameter :: MAXARGS = 50, MAXARGLEN = 30
	  integer :: nargsparsed
	  character(MAXARGLEN) :: arg(MAXARGS)

	  ! Integer list
	  integer, parameter :: MAXLISTITEMS = 20
	  type IntegerList
	    integer :: n = 0
	    integer, dimension(MAXLISTITEMS) :: items = 0
	  end type IntegerList
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

	! Parse string as comma-separated integer list into array (with index 0 containing the number of items)
	logical function parseIntegerList(string, list)
	implicit none
	character*(*) :: string
	type(IntegerList), intent(inout) :: list
	integer :: n, arglen, slen
	character :: i(MAXARGLEN)

	! Now get length of supplied string argument
	slen = len(string)
	!write(0,*) "String length = ", slen
	
	! Step through specified argument, looking for commas
	parseIntegerList = .true.
	arglen = 0
	i(:) = " "
	do n=1,slen
	  select case (ichar(string(n:n)))
	    case (44,32)
	      ! Comma (44) or space (32)
	      ! If arglen != 0 then 'store' current argument.
	      ! Otherwise, move on
	      if (arglen.ne.0) then
		list%n = list%n + 1
		if (list%n.gt.MAXLISTITEMS) then
		  write(0,*) "Integer list maximum size exceeded. MAXITEMS = ", MAXLISTITEMS
		  parseIntegerList = .false.
		  return
		end if
		read(i, "(i30)") list%items(list%n)
		arglen = 0
		i(:) = " "
		cycle
	      else if (ichar(string(n:n)).eq.44) then
		write(0,*) "!!! Found comma at start of argument?"
	        parseIntegerList = .false.
	        return
	      end if
	      exit
	    case (48:57)
	      ! Number - add to argument and continue
	      arglen = arglen + 1
	      i(arglen:arglen) = string(n:n)
	    case default
	      ! Not a comma, not a space, and not a number, so complain
	      write(0,*) "Found illegal character in integer list: ", string(n:n), " at pos ", n
	      parseIntegerList = .false.
	      return
	  end select
	end do

	end function parseIntegerList

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
