!	** moldist **
!	Calculate the distance between specified species(s), molecule(s), or atom(s)

	program moldist
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,a1,sp1,sp2,m2,mol2,baselen,nframes,success,nargs
	integer :: framestodo, count, start, finish
	integer :: iargc
	real*8 :: c1x,c1y,c1z,c2x,c2y,c2z,tx,ty,tz,rij,cutoff
	real*8, allocatable :: dist(:,:)

	nargs = iargc()
	if (nargs.LT.7) then
	  write(0,"(A)") "Usage : moldist <HISTORYfile> <OUTPUTfile> <sp1> <at1> <sp2> <mol2> <framestodo> [cutoff]"
	  write(0,"(A)") "        <sp1> <at1>      Central target - atom 'at1' in species 'sp1'"
	  write(0,"(A)") "        <sp2> <mol2>     Second species <sp2>. Specify specific molecule with <mol2>."
	  write(0,"(A)") "                            Otherwise, set <mol2> to 0 for all molecules of <sp2>."
	  write(0,"(A)") "        [cutoff]         If specified, do 'flagged' output for rij < cutoff."
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp); read(temp,"(I4)") sp1
	call getarg(4,temp); read(temp,"(I4)") a1
	call getarg(5,temp); read(temp,"(I4)") sp2
	call getarg(6,temp); read(temp,"(I4)") mol2
	call getarg(7,temp); read(temp,"(I6)") framestodo
	if (nargs.EQ.8) then
	  call getarg(8,temp); read(temp,"(F10.8)") cutoff
	else
	  cutoff = -1.0
	end if

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	allocate(dist(s_nmols(sp2),framestodo));
	dist(:,:) = reshape( (/ (0,n=1,s_nmols(sp2)*framestodo) /) , (/ s_nmols(sp2),framestodo /) )

	! XXXX
	! XXXX Main RDF routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes
	call calc_com

	! All distances will be calculated between the position of atom 'a1' in species 'sp1' and the com of 'sp2'
	c1x=xpos(s_start(sp1)-1+a1)
	c1y=ypos(s_start(sp1)-1+a1)
	c1z=zpos(s_start(sp1)-1+a1)
	finish = s_nmols(sp2); start = 1
	if (mol2.GT.0) then
	  start = mol2; finish = mol2
	end if
	count = 0
	do m2=start,finish
	  count = count + 1
	  ! Grab the centre-of-mass coordinates...
	  c2x=comx(sp2,m2)
	  c2y=comy(sp2,m2)
	  c2z=comz(sp2,m2)
	  ! Get the shortest (MIM) distance between the pair of positions...
	  call pbc(c2x,c2y,c2z,c1x,c1y,c1z,tx,ty,tz)
	  rij=sqrt( (tx-c1x)**2 + (ty-c1y)**2 + (tz-c1z)**2 )
	  ! Store this distance...
	  if (cutoff.GT.0) then
	    if (rij.LT.cutoff) dist(count,nframes) = m2
	  else
	    dist(count,nframes) = rij
	  end if
	end do

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

	hisfile = "moldist."
	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) then
	  basename="moldistresults."
	  baselen=11
	else
	  basename=hisfile(1:baselen)
	endif

	if (mol2.EQ.0) then
	  molpart = ""
	else
	  molpart = "m"//char(48+mol2/100)//char(48+(mol2-(mol2/100)*100)/10)//char(48+mod(mol2-(mol2/100)*100,10))
	end if
	if (cutoff.GT.0.0) then
	  resfile = basename(1:baselen)//"flag_s"//char(48+sp1)//"a"//char(48+a1/10)//char(48+MOD(a1,10))//"s"//char(48+sp2)//molpart
	else
	  resfile = basename(1:baselen)//"rij_s"//char(48+sp1)//"a"//char(48+a1/10)//char(48+MOD(a1,10))//"s"//char(48+sp2)//molpart
	end if
	open(unit=9,file=resfile,form="formatted",status="replace")
	finish = s_nmols(sp2); start = 1
	if (mol2.GT.0) then
	  start = mol2; finish = mol2
	end if
	count = 0
	do m2=start,finish
	  count = count + 1
	  do n=1,framestodo
	    if (cutoff.GT.0.0) then
	      write(9,"(4I6)") nint(dist(count,n))
	    else
	      write(9,"(F10.6)") dist(count,n)
	    end if
	  end do
	  write(9,*) ""
	  write(9,*) ""
	end do

	close(9)
	write(0,*) "Finished."
999	close(10)
	close(13)
	end program moldist
