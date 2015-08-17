	! Integer List

	module IList

	  ! Integer list
	  integer, parameter :: MAXLISTITEMS = 20
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
	integer :: n, arglen, slen
	character*20 :: i

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
	    case (43,45,48:57)
	      ! Plus, Minus, or Number - add to argument and continue
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
