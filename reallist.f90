	! Real List

	module RList

	  ! Real list
	  integer, parameter :: MAXLISTITEMS = 20
	  type RealList
	    integer :: n = 0
	    real, dimension(MAXLISTITEMS) :: items = 0
	  end type RealList

	contains

	! Parse string as comma-separated integer list
	logical function parseRealList(string, list)
	implicit none
	character*(*) :: string
	type(RealList), intent(inout) :: list
	integer :: n, arglen, slen
	character*20 :: i

	! Now get length of supplied string argument
	slen = len(string)
	!write(0,*) "String length = ", slen
	
	! Step through specified argument, looking for commas
	parseRealList = .true.
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
		  write(0,*) "Real list maximum size exceeded. MAXITEMS = ", MAXLISTITEMS
		  parseRealList = .false.
		  return
		end if
		read(i, "(i30)") list%items(list%n)
		arglen = 0
		i(:) = " "
		cycle
	      else if (ichar(string(n:n)).eq.44) then
		write(0,*) "!!! Found comma at start of argument?"
	        parseRealList = .false.
	        return
	      end if
	      exit
	    case (43,45,46,48:57,69,101)
	      ! Number - add to argument and continue
	      arglen = arglen + 1
	      i(arglen:arglen) = string(n:n)
	    case default
	      ! Not a comma, not a space, and not a number, so complain
	      write(0,*) "Found illegal character in integer list: ", string(n:n), " at pos ", n
	      parseRealList = .false.
	      return
	  end select
	end do

	end function parseRealList

	end module RList
