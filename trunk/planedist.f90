!	** planedist **
!	Calculate histogram of distances of atoms in species '2' from the z-plane defined for species 1

	program planedist
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,basename,resfile
	character*20 :: temp
	character*4 :: molpart
	integer :: n,m,o,s1,m1,m2,nframes,success,nargs,sp1,sp2, aoff, nbins
	integer :: framestodo = -1, framestoskip = 0, framesskipped, bin, baselen, nanglebins
	real*8, allocatable :: hist(:,:,:), angleacc(:,:)
	integer :: iargc
	real*8 :: tx,ty,tz,rij2,binwidth = 0.1, maxdist, rx, ry, rz, anglebin = 15.0
	real*8 :: v1(3), v2(3), xccangle, yccangle

	nargs = iargc()
	if (nargs.lt.4) then
	  write(0,"(A)") "Usage : planedist <DLP HISTORYfile> <DLP OUTPUTfile> <maxcomd> [-anglebin deg] [-skip n] [-frames n] [-axis sp x x y y] [-centresp sp] [-othersp sp] [-binwidth d]"
	  stop
	end if
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temp)
	read (temp,"(f10.5)") maxdist

	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
	allocate(aa(nspecies,4))

	n = 3
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-axis")
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") m
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(m,1)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(m,2)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(m,3)
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") aa(m,4)
	      write(0,"(A,I1,A,I2,A,I2,A,I2,A,I2,A)") "Molecular axis for species ",m," calculated from: X=",aa(m,1),"->", &
		 & aa(m,2),", Y=0.5X->0.5(r(",aa(m,3),")->r(",aa(m,4),"))"
	    case ("-centresp")
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp1
	      write(0,"(A,I4)") "Central species is = ",sp1
	    case ("-othersp")
	      n = n + 1; call getarg(n,temp); read(temp,"(I3)") sp2
	      write(0,"(A,I4)") "Other species is = ",sp2
	    case ("-frames")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestodo
	      write(0,"(A,I4)") "Frames to do is ", framestodo
	    case ("-skip")
	      n = n + 1; call getarg(n,temp); read(temp,"(I6)") framestoskip
	      write(0,"(A,I4)") "Frames to skip at start of trajectory is ", framestoskip
	    case ("-binwidth")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") binwidth
	      write(0,"(A,f10.5)") "Bin width to use is ", binwidth
	    case ("-anglebin")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") anglebin
	      write(0,"(A,f10.5)") "Y-COM-COM angle bin width to use is ", anglebin
	    case default
	      write(0,*) "Unrecognised argument: ", temp
	      stop
	  end select
	end do

	! Open and check the files...
	call openhis(hisfile,10)
	! Now, read in the history header so that we have cell()
	if (readheader().EQ.-1) goto 799

	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.EQ.-1) goto 799  ! File error....

	! Allocate arrays
	nbins = int(maxdist / binwidth) + 1
	nanglebins = int(180.0 / anglebin) + 1
	allocate(hist(s_natoms(sp2), 0:nanglebins, 0:nanglebins, 0:nbins))
	allocate(acc(0:nanglebins, 0:nanglebins))

	hist = 0.0
	acc = 0.0
	maxdist = maxdist * maxdist

	call alloc_axis

100	nframes=0
	framesskipped = 0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.EQ.-1) goto 799  ! File error....
	framesskipped = framesskipped + 1
	if (framesskipped.le.framestoskip) then
	  if (mod(framesskipped,100).EQ.0) write(0,*) "skipped... ",framesskipped
	  goto 101
	end if

	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	call calc_com
	call genaxis

	! For each atom of each molecule of species sp2, rotate into sp1's local frame and bin distance along z
	do m1=1,s_nmols(sp1)

	  aoff = s_start(sp2) - 1
	  do m2=1,s_nmols(sp2)

	    ! Check COM distance...
	    call pbc(comx(sp2,m2),comy(sp2,m2),comz(sp2,m2),comx(sp1,m1),comy(sp1,m1),comz(sp1,m1),tx,ty,tz)
	    tx = tx - comx(sp1,m1)
	    ty = ty - comy(sp1,m1)
	    tz = tz - comz(sp1,m1)
	    rij2 = tx*tx + ty*ty + tz*tz
	    if (rij2.le.maxdist) then

	      ! Now loop over atoms of species sp2 molecule m2
	      do n=1,s_natoms(sp2)

		! Get PBC atom delta with sp1/m1 COM
		call pbc(xpos(aoff+n),ypos(aoff+n),zpos(aoff+n),comx(sp1,m1),comy(sp1,m1),comz(sp1,m1),tx,ty,tz)
		tx = tx - comx(sp1,m1)
		ty = ty - comy(sp1,m1)
		tz = tz - comz(sp1,m1)
		rij2 = sqrt(tx*tx + ty*ty + tz*tz)

	write(0,*) rij2
		! Rotate this displacement t(xyz) into the local frame of sp1/m1
		rx = tx*axisx(sp1,m1,1) + ty*axisx(sp1,m1,2) + tz*axisx(sp1,m1,3)
		ry = tx*axisy(sp1,m1,1) + ty*axisy(sp1,m1,2) + tz*axisy(sp1,m1,3)
		rz = tx*axisz(sp1,m1,1) + ty*axisz(sp1,m1,2) + tz*axisz(sp1,m1,3)

		write(0,*) "rotated", sqrt(rx*rx+ry*ry+rz*rz)
	

		! Bin z-component
		bin = abs( int(rz / binwidth) )
		if (bin.le.nbins) then
		  hist(n,bin) = hist(n,bin) + 1.0
		end if

	      end do

	    end if

	    aoff = aoff + s_natoms(sp2)

	  end do !m2

	end do !m1

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

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  if (hisfile(n:n).EQ.".") THEN
	    baselen=n
	    goto 802
	  endif
	end do
802     if (baselen.EQ.-1) THEN
	  basename="rdfresults."
	  baselen=11
	ELSE
	  basename=hisfile(1:baselen)
	endif

        do n=1,s_natoms(sp2)
          resfile=basename(1:baselen)//"planedist"//CHAR(48+sp1)//"_"//CHAR(48+(n/10))//CHAR(48+MOD(n,10))
          OPEN(UNIT=9,file=resfile,FORM="FORMATTED")
          write(9,"(a,f12.4,i4,i4)") "#Bin     Hist    maxdist=",sqrt(maxdist),sp1,sp2
          ! Normalise histogram while writing...
          do m=0,nbins
            write(9,"(f6.3,3x,f12.8,4x,e12.6)") (m*binwidth)+binwidth/2.0, hist(n,m) / nframes
          end do
          close(9)
        end do

	write(0,*) "Finished."
999	close(10)
	close(13)
	end program planedist

