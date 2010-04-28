	! Parser Module

	module parse
	integer, parameter :: MAXARGS = 50
	integer :: nargsparsed
	character*30 :: arg(MAXARGS)

	contains

	! Read line from specified unit number and split into delimited arguments
	logical function readline(uno)
	implicit none
	character*200 line
	integer n,uno,arglen
	logical done
	! Blank previous arguments
	do n=1,MAXARGS
	  arg(n) = "                    "
	end do
	! Read line
	read(uno,"(a200)",end=900,err=900) line
	! Go through characters in line and get arguments
	nargsparsed = 0
	arglen = 0
	done = .false.
	do n=1,200
	  select case (line(n:n))
	    case (' ',',')
	      ! Blank character.
	      ! If arglen != 0 then 'store' current argument.
	      ! Otherwise, move on
	      if (arglen.ne.0) then
		arglen = arglen + 1
		! arg(nargsparsed+1)(arglen:arglen) = ENDOFLINE
		nargsparsed = nargsparsed + 1
		arglen = 0
	      end if
	    case ('!','#')
	      ! Comment character - we're done
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
	  read(arg(n),"(f20.10)") argr
	end if
	end function argr

	end module parse

