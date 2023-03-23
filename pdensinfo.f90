	! Display information on a pdens file

	program pdensinfo
	use parse; use PDensRW
	implicit none
	character*80 :: infile
	type(PDens) :: original
	integer :: iargc

        if (iargc().ne.1) stop "Usage: pdensinfo <pdensfile>"
        call getarg(1,infile)

	! Open pdens
	if (.not.loadPDens(infile, original)) stop "Failed to load pdens."

	end program pdensinfo
