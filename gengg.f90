	! Given a list of LJ params in a file, generates all cross-terms between unlike atoms
	! using geometric-geometric mixing rules

	program gengg
	use parse
	implicit none

	integer, parameter :: MAXTYPES=100
	character*80 :: paramfile
	real*8 :: sigma(MAXTYPES), eps(MAXTYPES)	! Storage for atomic sigma and epsilon from the file
	real*8 :: newsigma, neweps			! Calculated cross-term epsilon and sigma
	character*8 :: atomnames(MAXTYPES)		! Storage for atom type names from the file
	integer :: ntypes,n,m,i
	logical :: success
	integer :: iargc, nargs

	! Open file and read in data
	nargs = iargc()
	if (nargs.ne.1) stop "Usage:  gengg <paramfile>"
	call getarg(1,paramfile)
	open(unit=9, file='paramfile', form='formatted', status='old', err=999)

	ntypes = 0

10	success = readline(9)
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

	! Generate all cross terms using geometric combination rule
	! Determine total number of terms first
	i = 0
	do n=1,ntypes
	  i = i + n
	end do
	write(6,"('vdw ',i5)") i

	do n=1,ntypes
	  do m=n,ntypes
	    ! Write out data in DL_POLY format
	    newsigma = sqrt(sigma(n) * sigma(m))
	    neweps = sqrt(eps(n) * eps(m))
	    write(6,"(a8,a8,'lj',6x,f8.5,5x,f8.5)") atomnames(n), atomnames(m), neweps, newsigma
	  end do
	end do

	write(6,"('close',i5)") 

998	goto 1000
999	write(0,*) "Error."
	stop
1000	write(0,*) "Done."

	end program gengg
