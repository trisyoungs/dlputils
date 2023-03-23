!	** pairs **
!	Outputs XYZ coordinates of molecule pairs of two species that are within a supplied distance

	program pairs
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile,altheaderfile
	character*20 :: temp,xyzfile
	logical :: altheader = .FALSE., separate = .FALSE.
	integer :: n,m,s1,s2,m1,m2,baselen,nframes,success,o,nargs,numadded,i,j,framestodo = -1, u
	integer :: iargc, seed = -47849322
	real*8 :: dist, cutoff, tx, ty, tz, ranlimit = -1.0
	real*8, external :: ran2

	nargs = iargc()
	if (nargs.lt.5) stop "Usage : pairs <HISTORYfile> <OUTPUTfile> <sp1> <sp2> <dist> [-nframes n] [-header] [-ranlim r] [-seed i] [-separate]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
        call getarg(3,temp); read(temp,"(I6)") s1
        call getarg(4,temp); read(temp,"(I6)") s2
        call getarg(5,temp); read(temp,"(f10.6)") cutoff
	
	n = 5
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temp)
          select case (temp)
            case ("-nframes")
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
            case ("-separate")
              write(0,"(a,i4)") "Individual XYZ files will be written for each generated pair."
	      separate = .TRUE.
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
	resfile=basename(1:baselen)//"pairs"//CHAR(48+s1)//CHAR(48+s2)//".xyz"
	open(unit=15,file=resfile,form='formatted',status='replace')


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
	if (mod(nframes,100).EQ.0) write(0,"(i7,', numadded = ',i10)") nframes, numadded

	! Calculate molecule COMs
	call calc_com

	i = s_start(s1) - 1
	do m1=1,s_nmols(s1)
	  j = s_start(s2) - 1
	  do m2=1,s_nmols(s2)
	    ! Don't test distance if same species *and* same molecule
	    if ((m1.eq.m2).and.(s1.eq.s2)) cycle
	    ! Get the shortest (MIM) distance between the atom pair...
	    call pbc(comx(s2,m2),comy(s2,m2),comz(s2,m2),comx(s1,m1),comy(s1,m1),comz(s1,m1),tx,ty,tz)
	    dist=sqrt( (tx-comx(s1,m1))**2 + (ty-comy(s1,m1))**2 + (tz-comz(s1,m1))**2 )
	    ! Discard configurations randomly
	    if ((dist.lt.cutoff).and.(ran2(seed).gt.ranlimit)) then
	      numadded=numadded+1
	      ! Write out zero-centred and mim'd coordinates of the pair
	      ! For s1-m1 just subtract COM position from MIM atomic coordinates to this point
	      ! For s2-m2 get MIM coordinates w.r.t. COM of s1-m1, then subtract COM of s1-m1
	      do u=15,16
		if ((u.eq.16).and.(.not.separate)) cycle
		if (u.eq.16) then
		  xyzfile="pair     .xyz"
		  write(xyzfile(5:9),"(i5)") numadded
		  do n=5,9
		    if (xyzfile(n:n).eq." ") xyzfile(n:n) = "0"
		  end do
		  open(unit=16,file=xyzfile,form='formatted',status='replace')
		end if
	        write(u,"(i6)") s_natoms(s1)+s_natoms(s2)
	        write(u,"('Species ',i2,' Mol ',i4,' -- Species ',i2,' Mol ',i4,' -- Frame ',i4)") s1,m1,s2,m2,nframes
	        do n=1,s_natoms(s1)
		  call pbc(xpos(i+n),ypos(i+n),zpos(i+n),comx(s1,m1),comy(s1,m1),comz(s1,m1),tx,ty,tz)
		  write(u,"(a8,2x,3f12.6)") atmname(i+n)(1:1), tx-comx(s1,m1), ty-comy(s1,m1), tz-comz(s1,m1)
	        end do
	        do n=1,s_natoms(s2)
		  call pbc(xpos(j+n),ypos(j+n),zpos(j+n),comx(s1,m1),comy(s1,m1),comz(s1,m1),tx,ty,tz)
		  write(u,"(a8,2x,3f12.6)") atmname(j+n)(1:1), tx-comx(s1,m1), ty-comy(s1,m1), tz-comz(s1,m1)
	        end do
		if (u.eq.16) close(16)
	      end do
	    end if
	    j = j + s_natoms(s2)
	  end do
	  i = i + s_natoms(s1)
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
	end program pairs

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

