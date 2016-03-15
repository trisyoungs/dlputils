!	** pairs2 **
!	Outputs XYZ coordinates of first-shell molecule pairs of two species that are within a supplied distance

	program pairs2
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp,xyzfile
	logical :: altheader = .FALSE.
	integer :: n,m,sp1,sp2,sp,m1,m2,baselen,nframes,success,o,nargs,numadded,i,j,framestodo = -1, u
	integer :: iargc, seed = -47849322, setstodo, compairs(10,2), tth,th,hun,ten,units
	real*8 :: dist, mindist, maxdist, tx, ty, tz, ranlimit = -1.0, c1x, c1y, c1z, c2x, c2y, c2z
	real*8, external :: ran2

	nargs = iargc()
	if (nargs.lt.6) stop "Usage : pairs2 <HISTORYfile> <OUTPUTfile> <sp1> <sp2> <mindist> <maxdist> <nsets> [-frames n] [-header] [-ranlim r] [-seed i] [-compair sp i j]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
        call getarg(3,temp); read(temp,"(I6)") sp1
        call getarg(4,temp); read(temp,"(I6)") sp2
        call getarg(5,temp); read(temp,"(f10.6)") mindist
        call getarg(6,temp); read(temp,"(f10.6)") maxdist
        call getarg(7,temp); read(temp,"(I6)") setstodo

	compairs = 0

	n = 6
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-frames")
              n = n + 1; call getarg(n,temp); read(temp,"(i10)") framestodo
              write(0,"(a,i4)") "Frames to do:", framestodo
            case ("-seed")
              n = n + 1; call getarg(n,temp); read(temp,"(i10)") seed
              write(0,"(a,i10)") "Random seed is ",seed
            case ("-ranlim")
              n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") ranlimit
              write(0,"(a,f6.2)") "Random discard active - threshold = ",ranlimit
            case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(a,i4)") "Alternative header file supplied."
	      altheader = .TRUE.
            case ("-compair")
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") sp
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(sp,1)
              n = n + 1; call getarg(n,temp); read(temp,"(I6)") compairs(sp,2)
              write(0,"(A,3I4)") "Using COMpair for species ",sp, compairs(sp,:)
	  end select
	end do

	! Open and check the files...
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798

        ! Now, read in the history header so that we have cell()
	! If this fails then we may have a restarted trajectory. Continue, but only if
	! a header can be read in from the specified alternative history file..
	call openhis(hisfile,10)
        if (readheader().EQ.-1) then
	  if (altheader) then
	    write(0,*) "Restarted trajectory:"
	    close(dlpun_his)
	    call openhis(altheaderfile,10)
	    if (readheader().EQ.-1) goto 797
	    close(dlpun_his)
	    call openhis(hisfile,10)
	  else
	    goto 797
	  end if
	end if

	! Ascertain length of basename and open output file
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="pairs."
	  baselen=6
	else
	  basename=hisfile(1:baselen)
	endif

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
	numadded = 0
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,"(i7,', written ',i10,' of ',i10)") nframes, numadded, setstodo

	! Calculate molecule COMs
	call calc_com()

	! If compairs were specified, use that instead of COM
	do sp=1,nspecies
	  if (compairs(sp,1).ne.0) then
	    do m1=1,s_nmols(sp)
	      c1x = xpos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,1)-1)
	      c1y = ypos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,1)-1)
	      c1z = zpos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,1)-1)
	      c2x = xpos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,2)-1)
	      c2y = ypos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,2)-1)
	      c2z = zpos(s_start(sp)+(m1-1)*s_natoms(sp)+compairs(sp,2)-1)
	      call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	      comx(sp,m1) = (c1x+tx)*0.5
	      comy(sp,m1) = (c1y+ty)*0.5
	      comz(sp,m1) = (c1z+tz)*0.5
	    end do
	  end if
	end do

	i = s_start(sp1) - 1
	do m1=1,s_nmols(sp1)

	  ! Randomly select out central molecule
	  if (ran2(seed).lt.ranlimit) then
	    i = i + s_natoms(sp1)
	    cycle
	  end if

	  numadded=numadded+1

	  ! Find all close molecules of sp2 and write them out
	  m = 0
	  j = s_start(sp2) - 1
	  do m2=1,s_nmols(sp2)
	    ! Get the shortest (MIM) distance between the atom pair...
	    call pbc(comx(sp2,m2),comy(sp2,m2),comz(sp2,m2),comx(sp1,m1),comy(sp1,m1),comz(sp1,m1),tx,ty,tz)
	    dist=sqrt( (tx-comx(sp1,m1))**2 + (ty-comy(sp1,m1))**2 + (tz-comz(sp1,m1))**2 )

	    if ((dist.gt.mindist).and.(dist.lt.maxdist)) then

	      m = m + 1

	      ! Write out zero-centred and mim'd coordinates of all pairs with the distances specified...
	      ! For s1-m1 just subtract COM position from MIM atomic coordinates to this point
	      ! For s2-m2 get MIM coordinates w.r.t. COM of s1-m1, then subtract COM of s1-m1
	      ! Construct destfile name
	      tth = numadded / 10000; n = numadded - tth*10000
	      th = n / 1000; n = n - th*1000
	      hun = n / 100; n = n - hun*100
	      ten = n / 10; n = n - ten*10
	      units = n

	      resfile=basename(1:baselen)//"pairs"//".xyz"
	      resfile=basename(1:baselen)//"pairs"//CHAR(48+sp1)//CHAR(48+sp2)//"_set"//CHAR(48+tth)//CHAR(48+th)//CHAR(48+hun)//CHAR(48+ten)//CHAR(48+units)//"_pair"//CHAR(48+(m/10))//CHAR(48+MOD(m,10))//".xyz"
	      open(unit=12,file=resfile,form='formatted',status='replace')
	      write(12,"(i6)") s_natoms(sp1)+s_natoms(sp2)
	      write(12,"('Species ',i2,' Mol ',i4,' -- Species ',i2,' Mol ',i4,' -- Frame ',i4, 'dist=',f12.5)") sp1,m1,sp2,m2,nframes,dist
	      do n=1,s_natoms(sp1)
		call pbc(xpos(i+n),ypos(i+n),zpos(i+n),comx(sp1,m1),comy(sp1,m1),comz(sp1,m1),tx,ty,tz)
 		write(12,"(a8,2x,3f12.6)") atmname(i+n)(1:1), tx-comx(sp1,m1), ty-comy(sp1,m1), tz-comz(sp1,m1)
	      end do
	      do n=1,s_natoms(sp2)
		call pbc(xpos(j+n),ypos(j+n),zpos(j+n),comx(sp1,m1),comy(sp1,m1),comz(sp1,m1),tx,ty,tz)
		write(12,"(a8,2x,3f12.6)") atmname(j+n)(1:1), tx-comx(sp1,m1), ty-comy(sp1,m1), tz-comz(sp1,m1)
	      end do
	      close(12)
	    end if

	    j = j + s_natoms(sp2)
	  end do

	  if (numadded.eq.setstodo) goto 801
	  i = i + s_natoms(sp1)
	end do

	if (framestodo.eq.nframes) goto 801

	! Next frame
	goto 101

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

797	write(0,*) "No header found in history file. If a restarted trajectory, use '-header'"
	goto 999
798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended prematurely!"
	goto 801
800	write(0,*) "End of unformatted HISTORY file found."
801	write(0,*) ""

	write(0,*) "Finished."
999	close(15)
	end program pairs2

	real*8 function ran2(idum)
	implicit none
	integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
	real*8 :: AM,EPS,RNMX
	parameter (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1, &
	& IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
	& NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
	! 'ran2' random number generator.
	!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.
	integer :: idum2,j,k,iv(NTAB),iy
	save iv,iy,idum2
	data idum2/123456789/, iv/NTAB*0/, iy/0/
	if (idum.le.0) then
	  idum=max(-idum,1)
	  idum2=idum
	  do 11 j=NTAB+8,1,-1
	    k=idum/IQ1
	    idum=IA1*(idum-k*IQ1)-k*IR1
	    if (idum.lt.0) idum=idum+IM1
	    if (j.le.NTAB) iv(j)=idum
11	continue
	  iy=iv(1)
	endif
	k=idum/IQ1
	idum=IA1*(idum-k*IQ1)-k*IR1
	if (idum.lt.0) idum=idum+IM1
	k=idum2/IQ2
	idum2=IA2*(idum2-k*IQ2)-k*IR2
	if (idum2.lt.0) idum2=idum2+IM2
	j=1+iy/NDIV
	iy=iv(j)-idum2
	iv(j)=idum
	if(iy.lt.1)iy=iy+IMM1
	ran2=min(AM*iy,RNMX)
	return
	end function ran2

