	! Given a list of LJ params in a file, generates all cross-terms between unlike atoms

	program genlj
	use parse

	integer, parameter :: MAXTYPES=100
	real*8 :: sigma(MAXTYPES), eps(MAXTYPES)	! Storage for atomic sigma and epsilon from the file
	real*8 :: newsigma, neweps			! Calculated cross-term epsilon and sigma
	character*8 :: atomnames(MAXTYPES)		! Storage for atom type names from the file
	integer :: ntypes,n,m,i

	! Open file and read in data
	ntypes = 0

!10	read(5,"(a5,f10.5,f10.5)",end=99,err=99) atomnames(ntypes+1),eps(ntypes+1),sigma(ntypes+1)
10	success = readline(5)
	if (nargsparsed.le.0) goto 99
	ntypes = ntypes + 1
	atomnames(ntypes) = arg(1)(1:8)
	eps(ntypes) = argr(2)
	sigma(ntypes) = argr(3)
	goto 10

	! Print everything out to check
99	write(0,*) "ntypes ",ntypes
	write(0,*) "Type     Epsilon    Sigma"
	do n=1,ntypes
	  write(0,"(a8,2F10.5)") atomnames(n),eps(n),sigma(n)
	end do

	! Generate all cross terms using Lorentz-Berthelot combination rules
	! 'i' will count the number of terms we generate
	i = 0
	do n=1,ntypes
	  do m=n,ntypes
	    ! Write out data in DL_POLY format
	    newsigma = 0.5 * (sigma(n) + sigma(m))
	    neweps = sqrt(eps(n) * eps(m))
	    write(6,"(a8,a8,'lj',6x,f8.5,5x,f8.5)") atomnames(n), atomnames(m), neweps, newsigma
	    i = i + 1
	  end do
	end do

	write(0,*) "Terms generated:",i

	end program genlj
