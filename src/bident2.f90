!	** bident2 **
!	Calculate style of binding between 'bidentate' species and positions on a central species

	program bident2
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,headerfile
	character*20 :: temparg
	integer, parameter :: MAXSITES = 20, MAXCONTACTS = 10
	integer :: nframes,success,nargs,n,m1,m2,i,j,baselen,aoff1,aoff2,m
	integer :: framestodo = -1, framestodiscard = 0, sp1, sp2, nsp1sites = 0, sp1sites(MAXSITES,2), sp2sites(2)
	integer :: iargc, site, nsingle, nbridging, nbidentate, closebit, bit, isclose(MAXSITES)
	integer :: ncontacts
	integer, allocatable :: centralbits(:), sp2ncontacts(:), sp2contacts(:,:,:)
	logical :: dump = .false.
	real*8 :: a(3),b(3),c(3),mima(3),mimb(3),maxdist,dac,dbc, total, sitetotals(MAXSITES)
	real*8 :: singletab(MAXSITES), bifurtab(MAXSITES,MAXSITES), bridgetab(MAXSITES,MAXSITES), bitab(MAXSITES), multitab(MAXSITES)

	nargs = iargc()
	if (nargs.lt.7) then
	  write(0,"(a120)") "Usage : bident2 <HISTORYfile> <OUTfile> <sp1 (central)> <sp2 (outer)> <sp2 atom i> <sp2 atom j> <maxdist>"
	  write(0,"(a80)") "    [-frames n] [-discard n] [-site a b (e.g. H-O)] [-dump]"
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
              n = n + 1; call getarg(n,temparg); read(temparg,"(i6)") framestodiscard
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
	if (dump) then
	  open(unit=12,file=basename(1:baselen)//"bidentdump"//CHAR(48+sp1)//CHAR(48+sp2),form='formatted',status='replace')
	  write(12,*) "ncentre", s_nmols(sp1)
	  write(12,*) "nouter", s_nmols(sp2)
	  write(12,*) "ncentresites", nsp1sites
	end if

	allocate(centralbits(2**(nsp1sites*2)))
	allocate(sp2contacts(s_nmols(sp2),MAXCONTACTS,2))
	allocate(sp2ncontacts(s_nmols(sp2)))

	singletab = 0.0
	bitab = 0.0
	bridgetab = 0.0
	bifurtab = 0.0
	multitab = 0.0
	total = 0.0

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	if (framestodiscard.ne.0) then
	  framestodiscard = framestodiscard - 1
	  goto 101
	end if
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	if (dump) write(12,*) "frame", nframes

	sp2ncontacts = 0
	sp2contacts = 0

	! Inner loop over central species
	aoff1 = s_start(sp1)-1
	do m1 = 1,s_nmols(sp1)

	    ! Loop over outer species
	    aoff2 = s_start(sp2)-1
	    do m2=1,s_nmols(sp2)

	      a(1)=xpos(aoff2 + sp2sites(1))
	      a(2)=ypos(aoff2 + sp2sites(1))
	      a(3)=zpos(aoff2 + sp2sites(1))
	      b(1)=xpos(aoff2 + sp2sites(2))
	      b(2)=ypos(aoff2 + sp2sites(2))
	      b(3)=zpos(aoff2 + sp2sites(2))

	      isclose = 0
	      closebit = 0

	      ! Loop over sites on central species
	      do site = 1,nsp1sites

	        ! Grab this central site position
	        c(1)=xpos(aoff1+sp1sites(site,1))
	        c(2)=ypos(aoff1+sp1sites(site,1))
	        c(3)=zpos(aoff1+sp1sites(site,1))

	        ! Get MIM positions of outer species sites w.r.t. this site
	        call pbc(a(1),a(2),a(3),c(1),c(2),c(3),mima(1),mima(2),mima(3))
	        call pbc(b(1),b(2),b(3),c(1),c(2),c(3),mimb(1),mimb(2),mimb(3))
	        mima = mima - c
	        mimb = mimb - c

	        ! Determine distances
	        dac = sqrt(mima(1)*mima(1) + mima(2)*mima(2) + mima(3)*mima(3))
	        dbc = sqrt(mimb(1)*mimb(1) + mimb(2)*mimb(2) + mimb(3)*mimb(3))

	        ! Calculate closeness bit for central molecule
		bit = 0
	        if (dac.lt.maxdist) then
		  bit = bit + 2**(site-1)
		  isclose(site) = 1
		end if
		if (dbc.lt.maxdist) then
		  bit = bit + 2**(nsp1sites+site-1)
		  isclose(site) = isclose(site) + 2
		end if

	        if (dump.and.(bit.gt.0)) write(12,"(a,i5,a,i3,a,i5,a,i5)") "centre ",m1," site ", site, " contacts ",m2," bitmask ", isclose(site)
		closebit = closebit + bit

	      end do

	      ! Accumulate simple totals
	      j = 0
	      ncontacts = 0
	      do i=1,nsp1sites
	        if (isclose(i).eq.0) cycle
		if (isclose(i).lt.3) then
		  ncontacts = ncontacts + 1
		  j = i
		else
		  ncontacts = ncontacts + 2
		  j = i
		end if
	      end do

	      if (ncontacts.gt.0) then

		sp2ncontacts(m2) = sp2ncontacts(m2) + 1
		if (sp2ncontacts(m2).gt.MAXCONTACTS) stop "Exceeded MAXCONTACTS for outer species."
		sp2contacts(m2,sp2ncontacts(m2),1) = m1
		sp2contacts(m2,sp2ncontacts(m2),2) = closebit

	        select case (ncontacts)
		  case (1)
		    singletab(j) = singletab(j) + 1.0
		    !if (dump) write(12,*) "unique",j,m1
		  case (2)
		    ! Either a bidentate or bridging interaction
		    if (isclose(j).eq.3) then
		      bitab(j) = bitab(j) + 2.0
		      !if (dump) write(12,*) "bidentate",j,m1
		    else
		      ! Find other site...
		      do n=1,nsp1sites
		        if (n.eq.j) cycle
		        if (isclose(n).ne.0) exit
		      end do
		      ! Is it a bifurcated bridge or a 'bidentate' bridge?
		      if (isclose(n).eq.isclose(j)) then
		        bifurtab(n,j) = bifurtab(n,j) + 1.0
		        bifurtab(j,n) = bifurtab(j,n) + 1.0
		        !if (dump) write(12,*) "bifurcated",n,j,m1
		      else
		        bridgetab(n,j) = bridgetab(n,j) + 1.0
		        bridgetab(j,n) = bridgetab(j,n) + 1.0
		        !if (dump) write(12,*) "bridge",n,j,m1
		      end if
		    end if
		  case default
		    do i=1,nsp1sites
		      if (isclose(i).eq.0) cycle
		      if (isclose(i).lt.3) then
		        multitab(i) = multitab(i) + 1.0
		      else if (isclose(i).eq.3) then
		        multitab(i) = multitab(i) + 2.0
		      end if
		    end do
	        end select
	      end if

	      total = total + ncontacts
		    
	    aoff2 = aoff2 + s_natoms(sp2)

	  end do

	  aoff1 = aoff1 + s_natoms(sp1)

	end do

	if (dump) then
	  do m1=1,s_nmols(sp2)
	    if (sp2ncontacts(m1).eq.0) cycle
	    do n=1,sp2ncontacts(m1)
	      write(12,"(a,i5,a,i5,a,i5)") "outer ",m1," contacts ",sp2contacts(m1,n,1)," as ",sp2contacts(m1,n,2)
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

	! Determine site totals
	sitetotals = 0.0
	do i=1,nsp1sites
	  sitetotals(i) = sitetotals(i) + singletab(i)
	  sitetotals(i) = sitetotals(i) + bitab(i)
	  sitetotals(i) = sitetotals(i) + sum(bridgetab(i,:))
	  sitetotals(i) = sitetotals(i) + sum(bifurtab(i,:))
	  sitetotals(i) = sitetotals(i) + multitab(i)
	end do

	write(11,"('Found ',f12.2,' individual contacts in ', i7,' frames')") total, nframes
	write(11,"('  --> ',f12.2,' per frame')") total / nframes
	write(11,"('  --> ',f12.2,' per molecule')") total / nframes / s_nmols(sp1)
	write(11,"('Of these, ', f12.2, ' (',f5.1,'%) are isolated single contacts')") sum(singletab), sum(singletab)/total*100.0
	write(11,"('          --> ',f12.2,' per frame')") sum(singletab) / nframes
	write(11,"('          --> ',f12.2,' per central molecule')") sum(singletab) / nframes / s_nmols(sp1)
	write(11,"('          ', f12.2, ' (',f5.1,'%) are paired in bidentate (X-S-X) interactions')") sum(bitab), sum(bitab)/total*100.0
	write(11,"('          --> ',f12.2,' per frame')") sum(bitab) / nframes
	write(11,"('          --> ',f12.2,' per central molecule')") sum(bitab) / nframes / s_nmols(sp1)
	write(11,"('          ', f12.2, ' (',f5.1,'%) are pure bifurcated (S-X-S) interactions')") sum(bifurtab), sum(bifurtab)/total*100.0
	write(11,"('          --> ',f12.2,' per frame')") sum(bifurtab) / nframes
	write(11,"('          --> ',f12.2,' per central molecule')") sum(bifurtab) / nframes / s_nmols(sp1)
	write(11,"('          ', f12.2, ' (',f5.1,'%) are bridging (S-X-X-S) interactions')") sum(bridgetab), sum(bridgetab)/total*100.0
	write(11,"('          --> ',f12.2,' per frame')") sum(bridgetab) / nframes
	write(11,"('          --> ',f12.2,' per central molecule')") sum(multitab) / nframes / s_nmols(sp1)
	write(11,"('          ', f12.2, ' (',f5.1,'%) are multi-contact (3 or more) interactions')") sum(multitab), sum(multitab)/total*100.0
	write(11,"('          --> ',f12.2,' per frame')") sum(multitab) / nframes
	write(11,"('          --> ',f12.2,' per central molecule')") sum(multitab) / nframes / s_nmols(sp1)

	write(11,"(a)") ""
	write(11,"(a)") "Interaction breakdown for central molecule (per site per frame):"
	write(11,"(a)") ""
	write(11,"('Site       :', 20(3x,i2,4x))") (n,n=1,nsp1sites)
	write(11,"('Atom       :', 20(3x,i2,4x))") (sp1sites(n,1),n=1,nsp1sites)
	write(11,"(a)") ""
	write(11,"('N          :', 20(f7.3,2x))") (sitetotals(n)/nframes/s_nmols(sp1),n=1,nsp1sites)
	write(11,"('%All       :', 20(f7.3,2x))") (sitetotals(n)/total*100.0,n=1,nsp1sites)
	write(11,"(a)") ""
	write(11,"('Single     :', 20(f7.3,2x))") (singletab(n)/nframes/s_nmols(sp1),n=1,nsp1sites)
	write(11,"('%Group     :', 20(f7.3,2x))") (singletab(n)/sum(singletab)*100.0,n=1,nsp1sites)
	write(11,"('%Site      :', 20(f7.3,2x))") (singletab(n)/sitetotals(n)*100.0,n=1,nsp1sites)
	write(11,"('%All       :', 20(f7.3,2x))") (singletab(n)/total*100.0,n=1,nsp1sites)
	write(11,"(a)") ""
	write(11,"('Bidentate  :', 20(f7.3,2x))") (bitab(n)/nframes/s_nmols(sp1),n=1,nsp1sites)
	write(11,"('%Group     :', 20(f7.3,2x))") (bitab(n)/sum(bitab)*100.0,n=1,nsp1sites)
	write(11,"('%Site      :', 20(f7.3,2x))") (bitab(n)/sitetotals(n)*100.0,n=1,nsp1sites)
	write(11,"('%All       :', 20(f7.3,2x))") (bitab(n)/total*100.0,n=1,nsp1sites)
	write(11,"(a)") ""
	do m=1,nsp1sites
	  write(11,"('Bridge  ',i2,' :', 20(f7.3,2x))") m, (bridgetab(m,n)/nframes/s_nmols(sp1),n=1,nsp1sites)
	end do
	do m=1,nsp1sites
	  write(11,"('%Group  ',i2,' :', 20(f7.3,2x))") m, (bridgetab(m,n)/sum(bridgetab)*100.0,n=1,nsp1sites)
	end do
	write(11,"('%Site      :', 20(f7.3,2x))") (sum(bridgetab(n,:))/sitetotals(n)*100.0,n=1,nsp1sites)
	write(11,"('%All       :', 20(f7.3,2x))") (sum(bridgetab(n,:))/total*100.0,n=1,nsp1sites)
	write(11,"(a)") ""
	do m=1,nsp1sites
	  write(11,"('Bifurc  ',i2,' :', 20(f7.3,2x))") m, (bifurtab(m,n)/nframes/s_nmols(sp1),n=1,nsp1sites)
	end do
	do m=1,nsp1sites
	  write(11,"('%Group  ',i2,' :', 20(f7.3,2x))") m, (bifurtab(m,n)/sum(bifurtab)*100.0,n=1,nsp1sites)
	end do
	write(11,"('%Site      :', 20(f7.3,2x))") (sum(bifurtab(n,:))/sitetotals(n)*100.0,n=1,nsp1sites)
	write(11,"('%All       :', 20(f7.3,2x))") (sum(bifurtab(n,:))/total*100.0,n=1,nsp1sites)
	write(11,"(a)") ""
	write(11,"('Multi      :', 20(f7.3,2x))") (multitab(n)/nframes/s_nmols(sp1),n=1,nsp1sites)
	write(11,"('%Group     :', 20(f7.3,2x))") (multitab(n)/sum(multitab)*100.0,n=1,nsp1sites)
	write(11,"('%Site      :', 20(f7.3,2x))") (multitab(n)/sitetotals(n)*100.0,n=1,nsp1sites)
	write(11,"('%All       :', 20(f7.3,2x))") (multitab(n)/total*100.0,n=1,nsp1sites)

	close(11)

	write(0,*) "Finished."
999	close(14)
	end program bident2

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
        calcangle = safeAngle(dp)
	end function calcangle
