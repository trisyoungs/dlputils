!	** dahist
!	Calculate RDF between specific atoms of two species as well as angle map

	program dahist
	use dlprw; use utility
	implicit none
	integer, parameter :: maxtriplets = 50
	real*8, parameter :: pi = 3.14159265358979d0, radcon = 57.29577951d0
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	integer :: status,  nargs, success, baselen
	integer :: nbins,aoff1,aoff2,n,m,s1,s2,bin,angbin,nframes,numadded,sp1,sp2,a1(maxtriplets),a2(maxtriplets),a3(maxtriplets)
	integer :: iargc,p,framestodo,ntriplets
	real*8 :: r1(3), r2(3), r3(3), r12(3), r23(3), dist, binwidth, mindist, maxdist, angle, integral
	real*8, allocatable :: hist(:), anglehist(:), anglemap(:,:)

	binwidth=0.01   ! In Angstroms
	framestodo=-1
	sp1 = 1
	sp2 = 2
	mindist = 0.0
	ntriplets = 0

	nargs = iargc()
	if (nargs.lt.5) stop "Usage : dahist <HISTORYfile> <OUTPUTfile> <sp1> <sp2> <maxdist> [-triplet a1 a2 a3, where a1(sp1)...a2-a3(sp2)] [-frames n] [-mindist d]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read (temp,"(i4)") sp1
	call getarg(4,temp); read (temp,"(i4)") sp2
	call getarg(5,temp); read(temp,"(f20.14)") maxdist

	n = 5
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-frames") 
	      n = n + 1; call getarg(n,temp); read(temp,"(I5)") framestodo
	    case ("-mindist") 
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.5)") mindist
	    case ("-triplet") 
	      ntriplets = ntriplets + 1
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") a1(ntriplets)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") a2(ntriplets)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") a3(ntriplets)
	    case default
	      write(0,*) "Unrecognised command-line argument: ", temp
	      stop
	  end select
	end do
	     
	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	nbins = (maxdist-mindist) / binwidth
	write(0,"(A,i1,a,i1)") "'a1 / (a2,a3) taken from species ",sp1," / ",sp2
	write(0,"(A,F6.3,A)") "Using binwidth of ",binwidth," Angstroms"
	write(0,"(A,F6.3,A)") "Maximum a1...a2 distance to consider is ", maxdist
	write(0,"(A,F6.3,A)") "Minimum a1...a2 distance to consider is ", mindist
	write(0,"(A,I5,A)") "There will be ",nbins," histogram bins."
	if (framestodo.gt.0) write(0,"(a,i6)") "Number of frames to use in average = ",framestodo
	write(0,"(a)") "Triplets contributing to final maps are : "
	write(0,"(20(i3,'/ (',i3,',',i3,')'))") (a1(n),a2(n),a3(n),n=1,ntriplets)
	
	allocate(hist(nbins),stat=status); if (status.GT.0) stop "Allocation error for hist ()"
	allocate(anglehist(180),stat=status); if (status.GT.0) stop "Allocation error for anglehist ()"
	allocate(anglemap(nbins,180),stat=status); if (status.GT.0) stop "Allocation error for angle()"

	! Initialise the arrays...
	hist = 0.0
        anglehist = 0.0
	anglemap = 0.0

	! XXXX
	! XXXX Main loop
	! XXXX
	! Set up the vars...
	numadded = 0
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 120  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	
	aoff1 = s_start(sp1)
	do s1 = 1,s_nmols(sp1)     ! Loop over all molecules in species1

	  aoff2 = s_start(sp2) 
	  do s2 = 1,s_nmols(sp2)     ! Loop over all molecules of species 2

	    do n = 1,ntriplets

              ! Make sure to exclude i == k
              if ((aoff1+a1(n)).eq.(aoff2+a3(n))) cycle

	      ! Grab coordinates of what will be the central atom
	      r2(1) = xpos(aoff2+a2(n)-1)
	      r2(2) = ypos(aoff2+a2(n)-1)
	      r2(3) = zpos(aoff2+a2(n)-1)
 
	      ! Get minimum image positions of a1 and a3 w.r.t a2
	      call pbc(xpos(aoff1+a1(n)-1), ypos(aoff1+a1(n)-1), zpos(aoff1+a1(n)-1), r2(1), r2(2), r2(3), r1(1), r1(2), r1(3))
	      call pbc(xpos(aoff2+a3(n)-1), ypos(aoff2+a3(n)-1), zpos(aoff2+a3(n)-1), r2(1), r2(2), r2(3), r3(1), r3(2), r3(3))

	      r12 = r2 - r1
	      dist = sqrt( r12(1)*r12(1) + r12(2)*r12(2) + r12(3)*r12(3) )

	      if ((dist.ge.mindist).and.(dist.le.maxdist)) then
  
	        bin = int((dist-mindist) * (1.0/binwidth))+1
	        hist(bin) = hist(bin)+1.0

	        ! Determine angle
	        r23 = r2 - r3
	        r12 = r12 / dist
	        r23 = r23 / sqrt(sum(r23*r23))
	        angle = safeAngle(sum(r12*r23))
		write(0,*) aoff1+a1(n)-1, aoff2+a2(n)-1, aoff2+a3(n)-1, sum(r12*r23)

	        angbin = int(angle)+1
	        anglehist(angbin) = anglehist(angbin)+1.0

	        anglemap(bin,angbin) = anglemap(bin,angbin) + 1.0

	        numadded = numadded+1

	      end if

	    end do

	    aoff2 = aoff2 + s_natoms(sp2)
	  end do

	  aoff1 = aoff1 + s_natoms(sp1)

	end do   ! End main loop over all atoms of species1.

	if (nframes.EQ.1) then
	  write(0,"(A,I2,A,I4,A)") "RDF of atom in ",sp1," : averaged over ",s_nmols(sp1)," molecules."
	end if

	! Next frame
	if ((framestodo.gt.0).and.(framestodo.eq.nframes)) goto 120
	goto 101
	! Work on the results now to get the proper RDF
120	write(0,*) "Finished."
	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended prematurely!"
	goto 801
800	write(0,*) "End of unformatted HISTORY file found."

801	write(0,*) ""
	write(0,"(a,i6)") "Total frames used = ",nframes
	write(0,"(a,i10,a,i10)") "Total numadded = ",numadded,", per frame = ",numadded/nframes

	! Normalise data
	hist = hist / (nframes * s_nmols(sp1))
        anglehist = anglehist / sum(anglehist)
	anglemap = anglemap / (nframes * s_nmols(sp1))

	resfile=outputFileName(hisfile, "dahist", "dahist"//stringNMM(sp1,a1(1))//"_"//stringNMMOO(sp2,a2(1),a3(1))//".hist")
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,nbins
	  write(9,"(f6.3,3x,e16.8)") (n*binwidth+mindist)-binwidth/2.0,hist(n)
	end do
	close(9)

	resfile=outputFileName(hisfile, "dahist", "dahist"//stringNMM(sp1,a1(1))//"_"//stringNMMOO(sp2,a2(1),a3(1))//".angle")
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	do n=1,180
	  write(9,"(f6.1,3x,e16.8)") n-0.5,anglehist(n)
	end do
	close(9)

	resfile=outputFileName(hisfile, "dahist", "dahist"//stringNMM(sp1,a1(1))//"_"//stringNMMOO(sp2,a2(1),a3(1))//".surf")
	OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
	!write(9,*) nbins, 180
	!write(9,"(6f9.4)") binwidth,0.0,0.0,0.0,1.0,0.0
	!write(9,"(3f9.4)") 0.0, -1.0, 0.0
	!write(9,*) "yxz"
	do m=1,180
	  do n=1,nbins
	    write(9,"(3f16.8)") (n-0.5)*binwidth+mindist, (m-0.5)*1.0, anglemap(n,m)
	  end do
	  write(9,*) ""
	end do
	close(9)

	write(0,*) "Finished!"
999	CLOSE(10)
	CLOSE(13)
	end

