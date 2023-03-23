	! Integer List

	module IList

	  ! Integer list
	  integer, parameter :: MAXLISTITEMS = 50
	  type IntegerList
	    integer :: n = 0
	    integer, dimension(MAXLISTITEMS) :: items = 0
	  end type IntegerList

	contains

	! Parse string as comma-separated integer list
	logical function parseIntegerList(string, list)
	implicit none
	character*(*) :: string
	type(IntegerList), intent(inout) :: list
	integer :: n, arglen, slen, rangeStart, rangeEnd, rangeIndex
	logical :: range = .false.
	character*20 :: i, j

	! Now get length of supplied string argument
	slen = len(string)
	
	! Step through specified argument, looking for commas
	parseIntegerList = .true.
	arglen = 0
	i(:) = " "
	j(:) = " "
	do n=1,slen
	  select case (ichar(string(n:n)))
	    case (44,32)
	      ! Comma (44) or space (32)
	      ! If arglen != 0 then 'store' current argument.
	      ! Otherwise, move on
	      if (arglen.ne.0) then
		! If its a range, add each integer in turn
		if (range) then
		  read(i, "(i30)") rangeEnd
		  read(j, "(i30)") rangeStart
		  do rangeIndex=rangeStart,rangeEnd
		    list%n = list%n + 1
		    if (list%n.gt.MAXLISTITEMS) then
		      write(0,*) "Integer list maximum size exceeded. MAXITEMS = ", MAXLISTITEMS
		      parseIntegerList = .false.
		      return
		    end if
		    list%items(list%n) = rangeIndex
		  end do
		else
		  list%n = list%n + 1
		  if (list%n.gt.MAXLISTITEMS) then
		    write(0,*) "Integer list maximum size exceeded. MAXITEMS = ", MAXLISTITEMS
		    parseIntegerList = .false.
		    return
		  end if
		  read(i, "(i30)") list%items(list%n)
		end if
		arglen = 0
		i(:) = " "
		j(:) = " "
		range = .false.
		cycle
	      else if (ichar(string(n:n)).eq.44) then
		write(0,*) "!!! Found comma at start of argument?"
	        parseIntegerList = .false.
	        return
	      end if
	      exit
	    case (45)
	      ! Minus - could be a range specified, so check for zero argument length (in which case it's a unary minus)
	      if (arglen.eq.0) then
	        arglen = arglen + 1
	        i(arglen:arglen) = '-'
	      else
		! It's a range, so store current argument i in j, and reset to continue
		j = i
		i(:) = " "
		arglen = 0
		range = .true.
	      end if
	    case (43,48:57)
	      ! Plus, or Number - add to argument and continue
	      arglen = arglen + 1
	      i(arglen:arglen) = string(n:n)
	    case (88,120)
	      ! Upper or lowercase 'x' - terminate parsing
	      return
	    case default
	      ! Not a comma, not a space, and not a number, so complain
	      write(0,*) "Found illegal character in integer list: ", string(n:n), " at pos ", n
	      parseIntegerList = .false.
	      return
	  end select
	end do

	end function parseIntegerList

	! Return whether specified value is contained within the list
	logical function integerListContains(list, item)
	implicit none
	type(IntegerList), intent(in) :: list
	integer, intent(in) :: item
	integer :: n
	
	integerListContains = .false.
	do n=1,list%n
	  if (list%items(n).eq.item) then
	    integerListContains = .true.
	    return
	  end if
	end do

	end function integerListContains

	end module IList
