	! Calculate statistics for supplied data
	program stats
	use parse
	implicit none
	integer :: ndat, n, success
	real*8 :: dat(100000), avg, sd, minx, maxx, sumxsq

	! Read in data
	ndat = 0
10	success = readline(5)
	if (nargsparsed.eq.-1) goto 100
	if (nargsparsed.eq.0) goto 10
	ndat = ndat + 1
	dat(ndat) = argr(1)
	goto 10

100	write(6,*) "Read ",ndat,"data values"

        ! Calculate range values
	minx = dat(1)
	maxx = dat(1)
	avg = 0.0d0
        do n=1,ndat
          minx = min(minx,dat(n))
          maxx = max(maxx,dat(n))
          avg = avg + dat(n)
	end do
	avg = avg / real(ndat)

	! Calculate standard deviation
        sumxsq = 0.0d0
	do n=1,ndat
          sumxsq = sumxsq + (dat(n) - avg)**2
	end do
        sd = SQRT( sumxsq / real(ndat) )

	! Write results (avg, sd, min, max)
	write(6,"(a)") "     Average         S.D.        Minimum       Maximum"
	write(6,"(4e14.6)") avg, sd, minx, maxx

	end program stats
