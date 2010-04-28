!	** bident **
!	Calculate style of binding between 'bidentate' species and positions on a central species

	program bident
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,headerfile
	character*20 :: temparg
	integer, parameter :: MAXSITES = 20
	integer :: nframes,success,nargs,n,m1,m2,i,j,baselen,aoff1,aoff2,m
	integer :: framestodo = -1, framestoskip = 0, sp1, sp2, nsp1sites = 0, sp1sites(MAXSITES,2), sp2sites(2)
	integer :: iargc, site
	logical :: aclose(MAXSITES), bclose(MAXSITES), dump = .false.
	real*8 :: a(3),b(3),c(3),mima(3),mimb(3),maxdist,dac,dbc, closetab(MAXSITES), bridgetab(MAXSITES,MAXSITES), bitab(MAXSITES)

	nargs = iargc()
	if (nargs.lt.7) then
	  write(0,"(a120)") "Usage : bident <HISfile> <OUTfile> <sp1 (central)> <sp2 (outer)> <sp2 atom i> <sp2 atom j> <maxdist>"
	  write(0,"(a80)") "    [-nframes n] [-discard n] [-site a b (e.g. H-O)] [-dump]"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temparg); read(temparg,"(I5)") sp1
	call getarg(4,temparg); read(temparg,"(I5)") sp2
	call getarg(5,temparg); read(temparg,"(I5)") sp2sites(1)
	call getarg(6,temparg); read(temparg,"(I5)") sp2sites(2)
	call getarg(7,temparg); read(temparg,"(f10.5)") maxdist
	
	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	n = 7
        do
          n = n + 1; if (n.GT.nargs) exit
          call getarg(n,temparg)
          select case (temparg)
            case ("-site")
	      nsp1sites = nsp1sites + 1
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") sp1sites(nsp1sites,1)
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") sp1sites(nsp1sites,2)
              write(0,"(A,i4,' and ',i4)") "Site on species 1 defined: atoms are ",sp1sites(nsp1sites,1),sp1sites(nsp1sites,2)
            case ("-header")
              n = n + 1; call getarg(n,headerfile)
            case ("-frames")
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") framestodo
            case ("-discard")
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") framestoskip
            case ("-dump")
              dump = .true.
	    case default
	      write(0,*) "Unrecognised argument:", temparg
	      stop
	  end select
	end do
	
	! Open output files
	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="rdfresults."
	  baselen=11
	else
	  basename=hisfile(1:baselen)
	endif
	open(unit=11,file=basename(1:baselen)//"bident"//CHAR(48+sp1)//CHAR(48+sp2),form='formatted',status='replace')
	write(11,"('# Central species is ',i4)") sp1
	write(11,"('# Number of sites defined on central species ',i4)") nsp1sites
	do n=1,nsp1sites
	  write(11,"(6x,2i5)") sp1sites(n,1),sp1sites(n,2)
	end do
	write(11,"('# Outer species is ',i4)") sp2
	write(11,"('# Atom IDs used for interaction points :',i4,i4)") sp2sites
	if (dump) open(unit=12,file=basename(1:baselen)//"dump"//CHAR(48+sp1)//CHAR(48+sp2),form='formatted',status='replace')

	closetab = 0.0
	bitab = 0.0
	bridgetab = 0.0

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.EQ.-1) goto 799  ! File error....
	if (framestoskip.ne.0) then
	  framestoskip = framestoskip - 1
	  goto 101
	end if
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	if (dump) write(12,*) "frame", nframes

	! Outer loop over 'outer' species
	aoff2 = s_start(sp2)-1
	do m2=1,s_nmols(sp2)

	  ! Grab coordinates of the outer species two sites
	  a(1)=xpos(aoff2 + sp2sites(1))
	  a(2)=ypos(aoff2 + sp2sites(1))
	  a(3)=zpos(aoff2 + sp2sites(1))
	  b(1)=xpos(aoff2 + sp2sites(2))
	  b(2)=ypos(aoff2 + sp2sites(2))
	  b(3)=zpos(aoff2 + sp2sites(2))
	  
	  ! Inner loop over central species
	  aoff1 = s_start(sp1)-1
	  do m1 = 1,s_nmols(sp1)

	    ! Loop over sites on central species
	    do site = 1,nsp1sites

	      ! Grab this central site position
	      c(1)=xpos(aoff1+sp1sites(site,1))
	      c(2)=ypos(aoff1+sp1sites(site,1))
	      c(3)=zpos(aoff1+sp1sites(site,1))

	      ! Get MIM positions of sp2 sites w.r.t. this site
	      call pbc(a(1),a(2),a(3),c(1),c(2),c(3),mima(1),mima(2),mima(3))
	      call pbc(b(1),b(2),b(3),c(1),c(2),c(3),mimb(1),mimb(2),mimb(3))
	      mima = mima - c
	      mimb = mimb - c

	      ! Determine distances
	      dac = sqrt(mima(1)*mima(1) + mima(2)*mima(2) + mima(3)*mima(3))
	      dbc = sqrt(mimb(1)*mimb(1) + mimb(2)*mimb(2) + mimb(3)*mimb(3))

	      ! Calculate closeness
	      aclose(site) = dac.lt.maxdist
	      bclose(site) = dbc.lt.maxdist

	    end do

	    ! Analyse data
	    do i=1,nsp1sites
	      if (.not.aclose(i)) cycle
	      closetab(i) = closetab(i) + 1.0
	      do j=1,nsp1sites
		if (.not.bclose(j)) cycle
	        closetab(j) = closetab(j) + 1.0
		if (i.ne.j) bridgetab(i,j) = bridgetab(i,j) + 1.0
		if (i.ne.j) bridgetab(j,i) = bridgetab(j,i) + 1.0
		if (i.eq.j) bitab(i) = bitab(i) + 1.0
	      end do
	    end do

	    aoff1 = aoff1 + s_natoms(sp1)

	  end do

	  aoff2 = aoff2 + s_natoms(sp2)

	end do

	! Write dump file info
	if (dump) then
	  do m1=1,s_nmols(sp1)
	    write(12,*) "mol", m1
	    ! For each central molecule, determine which interactions exist....
	    do site=1,nsp1sites
	    end do
	  end do
	end if

	if (nframes.EQ.framestodo) goto 800
	! Next frame
	goto 101

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999
798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended before framestodo was fulfilled..."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframes," (wanted ",framestodo,")"
	goto 801
800	write(0,*) "Framestodo was fulfilled."
801	write(0,*) ""

	write(11,"('Found ',f12.1,' individual interactions in ', i7,' frames')") sum(closetab), nframes
	write(11,"('  --> ',f12.1,' per frame')") sum(closetab) / nframes
	write(11,"('  --> ',f12.1,' per molecule')") sum(closetab) / nframes / s_nmols(sp1)
	write(11,"('Of these, ', f12.1, ' (',f5.1,'%) are accounted for by bidentate interactions')") sum(bitab)*2.0, sum(bitab)*2.0/sum(closetab)*100.0
	write(11,"('  --> ',f12.1,' per frame')") sum(bitab)*2.0 / nframes
	write(11,"('  --> ',f12.1,' per molecule')") sum(bitab)*2.0 / nframes / s_nmols(sp1)

	write(11,"(a)") ""
	write(11,"(a)") "Interaction counts (per molecule per frame):"
	write(11,"(a)") ""
	write(11,"('Site    :', 20(3x,i2,4x))") (n,n=1,nsp1sites)
	write(11,"('Atom    :', 20(3x,i2,4x))") (sp1sites(n,1),n=1,nsp1sites)
	write(11,"(a)") ""
	write(11,"('N       :', 20(f7.3,2x))") (closetab(n)/nframes/s_nmols(sp1),n=1,nsp1sites)
	write(11,"('%All    :', 20(f7.3,2x))") (closetab(n)/sum(closetab)*100.0,n=1,nsp1sites)
	write(11,"(a)") ""
	write(11,"('Single  :', 20(f7.3,2x))") ((closetab(n)-bitab(n)*2-sum(bridgetab(n,:)))/nframes/s_nmols(sp1),n=1,nsp1sites)
	write(11,"('%Site   :', 20(f7.3,2x))") ((closetab(n)-bitab(n)*2-sum(bridgetab(n,:)))/closetab(n)*100.0,n=1,nsp1sites)
	write(11,"(a)") ""
	write(11,"('Bi      :', 20(f7.3,2x))") (bitab(n)*2/nframes/s_nmols(sp1),n=1,nsp1sites)
	write(11,"('%Bi     :', 20(f7.3,2x))") (bitab(n)/sum(bitab)*100.0,n=1,nsp1sites)
	write(11,"('%Site   :', 20(f7.3,2x))") (bitab(n)*2/closetab(n)*100.0,n=1,nsp1sites)
	write(11,"(a)") ""
	do m=1,nsp1sites
	  write(11,"('Br   ',i2,' :', 20(f7.3,2x))") m, (bridgetab(m,n)/nframes/s_nmols(sp1),n=1,nsp1sites)
	end do
	do m=1,nsp1sites
	  write(11,"('%Br  ',i2,' :', 20(f7.3,2x))") m, (bridgetab(m,n)/sum(bridgetab)*100.0,n=1,nsp1sites)
	end do
	write(11,"('%Site   :', 20(f7.3,2x))") (sum(bridgetab(n,:))/closetab(n)*100.0,n=1,nsp1sites)

	close(11)

	write(0,*) "Finished."
999	close(14)
	end program bident

	real*8 function calcangle(a,b,c)
	use utility; implicit none
	real*8,dimension(3) :: a,b,c,v1,v2,temp
	real*8 :: v1mag,v2mag,dp
        ! Angle between atoms a--b--c
	call pbc(a(1),a(2),a(3),b(1),b(2),b(3),temp(1),temp(2),temp(3))
	a = temp
	call pbc(c(1),c(2),c(3),b(1),b(2),b(3),temp(1),temp(2),temp(3))
	c = temp
	! Get the bond vectors
	v1 = b - a
	v2 = b - c
	v1mag = sqrt( v1(1)**2 + v1(2)**2 + v1(3)**2 )
	v2mag = sqrt( v2(1)**2 + v2(2)**2 + v2(3)**2 )
	v1 = v1 / v1mag
	v2 = v2 / v2mag
        ! Calculate dot product and angle...
        dp = (v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3))
        calcangle = acos(dp)* 57.29577951d0
	end function calcangle
