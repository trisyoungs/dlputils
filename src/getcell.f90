!	** getcell **
!	Retrieves the cell dimensions from a HISTORY file and outputs them to stdout

	program getcell
	use dlprw; use utility; implicit none
	character*80 :: hisfile, headerfile, outfile
	character*8 :: discard
	integer :: success, nargs
	integer :: iatm, n, discardn, nframes, l, fileform
	integer :: iargc
	logical :: altheader = .false.
	real*8, parameter :: radcon = 57.29577951d0
	! Data is lengths(3) - angles(3) - volume
	real*8 :: celldata(7), mindata(7), maxdata(7), avgdata(7)

	nargs = iargc()
	if (nargs.lt.1) stop "Usage : getcell <HISTORYfile> [HISTORY altheader]"
        call getarg(1, hisfile)
	if (nargs.eq.2) then
	  altheader = .true.
          call getarg(2, headerfile)
	end if

	write(6,"(A)") "#    Step           A               B             C         Alpha      Beta     Gamma    Volume   "

	! Set the necessary variables to zero.....
	mindata = 1e8
	maxdata = -1e8
	avgdata = 0.0
	nframes=0

        ! Read the header of the history file...
        call openhis(hisfile,15)
        if (readheader().EQ.-1) then
          if (altheader) then
            write(0,*) "Restarted trajectory:"
            close(dlpun_his)
            call openhis(headerfile,10)
            if (readheader().EQ.-1) stop "Couldn't read header of alternative file."
            close(dlpun_his)
            call openhis(hisfile,10)
          else
            stop "Couldn't read header of history file."
          end if
        end if

	! At START of trajectory
100	success=readframe()
	if (success.EQ.1) goto 800  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	!if (MOD(nframes,100).EQ.0) write(0,*) nframes

	! Calculate cell lengths...
	celldata(1) = sqrt( cell(1)*cell(1) + cell(2)*cell(2) + cell(3)*cell(3) )
	celldata(2) = sqrt( cell(4)*cell(4) + cell(5)*cell(5) + cell(6)*cell(6) )
	celldata(3) = sqrt( cell(7)*cell(7) + cell(8)*cell(8) + cell(9)*cell(9) )
	! ...cell volume (matrix determinant)...
	celldata(7) = cell(1) * (cell(5)*cell(9) - cell(8)*cell(6));
	celldata(7) = celldata(7) - cell(2) * (cell(4)*cell(9) - cell(7)*cell(6));
	celldata(7) = celldata(7) + cell(3) * (cell(4)*cell(8) - cell(7)*cell(5));
	! ...and cell angles
	do n=1,3
	  cell(n) = cell(n) / celldata(1)
	  cell(n+3) = cell(n+3) / celldata(2)
	  cell(n+6) = cell(n+6) / celldata(3)
	end do
	celldata(4) = cell(4)*cell(7) + cell(5)*cell(8) + cell(6)*cell(9)
	celldata(5) = cell(1)*cell(7) + cell(2)*cell(8) + cell(3)*cell(9)
	celldata(6) = cell(1)*cell(4) + cell(2)*cell(5) + cell(3)*cell(6)
	celldata(4) = safeAngle(celldata(4))
	celldata(5) = safeAngle(celldata(5))
	celldata(6) = safeAngle(celldata(6))

	! Store min, max, and average
	do n=1,7
	  if (celldata(n).lt.mindata(n)) mindata(n) = celldata(n)
	  if (celldata(n).gt.maxdata(n)) maxdata(n) = celldata(n)
	  avgdata(n) = avgdata(n) + celldata(n)
	end do

	! Now write out the data.....
	write(6,"(i8,6(2x,f13.8),2x,f13.4)") nstep,celldata

	goto 100

799	write(0,*) "HISTORY file ended prematurely!"
	goto 805
800	write(6,*) "End of unformatted HISTORY file found."
805	close(15)

	write(0,*) "Final averages......"
	avgdata = avgdata / nframes
	write(0,"(a8,6(2x,f13.8),2x,f13.4)") "Minimum",mindata
	write(0,"(a8,6(2x,f13.8),2x,f13.4)") "Average",avgdata
	write(0,"(a8,6(2x,f13.8),2x,f13.4)") "Maximum",maxdata
	write(6,"('#',a7,6(2x,f13.8),2x,f13.4)") "Minimum",mindata
	write(6,"('#',a7,6(2x,f13.8),2x,f13.4)") "Average",avgdata
	write(6,"('#',a7,6(2x,f13.8),2x,f13.4)") "Maximum",maxdata
	
999	close(15)
	end program getcell

