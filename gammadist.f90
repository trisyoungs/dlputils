!	** gammadist.f90 **
!	Program to calculate the vector of molecules (pointing vector) in a system

	program gammadist
	use dlprw; use utility
	implicit none
	integer, parameter :: NGRID=100
	integer :: maxframes,baselen,nargs,success,n,i,targetsp,startf,endf,oid,hid(2)
	integer :: nframes, nframesused, xbin, ybin, zspecies, surfsp, surfid
	character*80 :: hisfile,outfile,basename,resfile,temp,flagfile,altheaderfile
	character*8 :: discard
	logical :: altheader = .FALSE.
	real*8 :: acc(-NGRID:NGRID,-NGRID:NGRID) 
	real*8 :: px,py,pz, tx,ty,tz, dist, 
	real*8 :: map(-NGRID:NGRID,-NGRID:NGRID,5)
	real*8 :: minz, maxz, xdelta, ydelta, zorigin
	integer :: iargc

	write(0,*) "*** pdens"

	nargs = iargc()
	if (nargs.lt.8) then
	  write(*,"(a)") "Usage: gammadist <HISTORYfile> <OUTPUTfile> <watersp> <O id> <H1 id> <H2 id> <surfsp> <surfid> [options]"
	  write(*,"(a)") "        [-zorigin z]            Use specified Z as first distance point"
	  write(*,"(a)") "        [-minz r]               Minimum deviation from first distance point"
	  write(*,"(a)") "        [-maxz r]               Maximum deviation from first distance point"
	  write(*,"(a)") "        [-header file]          Use specified file to get header"
	  stop
	else
	  call getarg(1,hisfile)
	  call getarg(2,outfile)
	  call getarg(3,temp); read(temp,"(i4)") targetsp
	  call getarg(4,temp); read(temp,"(i4)") oid
	  call getarg(5,temp); read(temp,"(i4)") hid(1)
	  call getarg(6,temp); read(temp,"(i4)") hid(2)
	  call getarg(7,temp); read(temp,"(i4)") surfsp
	  call getarg(8,temp); read(temp,"(i4)") surfid
	end if

	! Ascertain length of basename....
	baselen=-1
	do n=80,1,-1
	  IF (hisfile(n:n).eq.".") then
	    baselen=n
	    goto 50
	  endIF
	end do
50	if (baselen.eq.-1) then
	  basename="pointresults."
	  baselen=14
	else
	  basename=hisfile(1:baselen)
	endif

	open(unit=20,file=basename(1:baselen)//"out",status='replace',form='formatted')

	write(0,"(A,A)") "History file : ",hisfile
	write(20,"(A,A)") "History file : ",hisfile
	!call openhis(hisfile,n)
	write(0,"(A,A)") " Output file : ",outfile
	write(20,"(A,A)") " Output file : ",outfile
	if (outinfo(outfile,1).eq.-1) goto 798

	! Set some variable defaults before we read in any command-line arguments
	minz = 0.0
	maxz = 100.0
	zorigin = 0.0
	map = 0.0
	acc = 0.0

	n = 2
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temp)
	  select case (temp)
	    case ("-maxz")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") maxz
	      write(0,"(A,f8.4)") "Maximum z deviation = ",maxz
	      write(20,"(A,f8.4)") "Maximum z deviation = ",maxz
	    case ("-minz")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") minz
	      write(0,"(A,f8.4)") "Minimum z deviation = ",minz
	      write(20,"(A,f8.4)") "Minimum z deviation = ",minz
	    case ("-header")
              n = n + 1; call getarg(n,altheaderfile)
              write(0,"(A)") "Alternative header file supplied."
              write(20,"(A)") "Alternative header file supplied."
              altheader = .TRUE.
	    case ("-zorigin")
	      n = n + 1; call getarg(n,temp); read(temp,"(f10.4)") zorigin
	      write(0,"(A,f8.4)") "Z-basis origin is = ",zorigin
	      write(20,"(A,f8.4)") "Z-basis origin is = ",zorigin
	    case default
	      write(0,*) "Unrecognised command line option : ",temp
	      stop
	  end select
	end do
	
	! Read the header of the history file...
        call openhis(hisfile,10)
        if (readheader().EQ.-1) then
          if (altheader) then
            write(0,*) "Restarted trajectory:"
            close(dlpun_his)
            call openhis(altheaderfile,10)
            if (readheader().EQ.-1) goto 799
            close(dlpun_his)
            call openhis(hisfile,10)
          else
            goto 799
          end if
        end if

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes = 0
	nframesused = 0
	xdelta = cell(1) / real(2*NGRID+1)
	ydelta = cell(5) / real(2*NGRID+1)

102	success=readframe()
	IF (success.eq.1) goto 801  ! End of file encountered....
	IF (success.eq.-1) goto 119  ! File error....
	nframes=nframes+1
	if (MOD(nframes,100).eq.0) write(0,*) nframes
	if (nframes.LT.startf) goto 102

	nframesused = nframesused + 1

	! Loop over 'water' molecules
	i = s_start(targetsp)
	do n=1,s_nmols(targetsp)

	  ! Loop over H atoms
	  do h=1,2

	    ! Check this H atom is in the specified z-slice
	    if ((dabs(tz - zorigin).gt.maxz).or.(dabs(tz-zorigin).lt.minz)) then
	      i = i + s_natoms(targetsp)
	      cycle
	    end if



	  end do



	  ! Get normalised MIM vector of first atom with third atom (in tx)
	  call pbc(xpos(i+vec(2)-1),ypos(i+vec(2)-1),zpos(i+vec(2)-1), &
		& xpos(i+vec(3)-1),ypos(i+vec(3)-1),zpos(i+vec(3)-1), tx, ty, tz)
	  mag = sqrt(tx*tx + ty*ty + tz*tz)
	  tx = tx / mag
	  ty = ty / mag
	  tz = tz / mag
	  
	  ! Get MIM vector of second atom with third atom (in px)
	  call pbc(xpos(i+vec(2)-1),ypos(i+vec(2)-1),zpos(i+vec(2)-1), &
		& xpos(i+vec(3)-1),ypos(i+vec(3)-1),zpos(i+vec(3)-1), px, py, pz)
	  mag = sqrt(px*px + py*py + pz*pz)
	  px = px / mag
	  py = py / mag
	  pz = pz / mag

	  ! Get average of these two vectors (in tx)
	  tx = (tx + px) * 0.5
	  ty = (ty + py) * 0.5
	  tz = (tz + pz) * 0.5

	  ! Calculate and normalise pointing vector
	  pvec(1) = xpos(i+vec(3)-1) - tx
	  pvec(2) = ypos(i+vec(3)-1) - ty
	  pvec(3) = zpos(i+vec(3)-1) - tz
	  mag = sqrt(pvec(1)*pvec(1) + pvec(2)*pvec(2) + pvec(3)*pvec(3))
	  pvec(1) = pvec(1) / mag
	  pvec(2) = pvec(2) / mag
	  pvec(3) = pvec(3) / mag

	  ! Get x/y position on surface
	  ! Fold into -0.5->0.5
	  if (tx.lt.(-0.5*cell(1))) tx = tx + cell(1)
	  if (tx.gt.(0.5*cell(1))) tx = tx - cell(1)
	  if (ty.lt.(-0.5*cell(5))) ty = ty + cell(5)
	  if (ty.gt.(0.5*cell(5))) ty = ty - cell(5)
	  if (tz.lt.(-0.5*cell(9))) tz = tz + cell(9)
	  if (tz.gt.(0.5*cell(9))) tz = tz - cell(9)
	  ! Check z-distance
	  if ((dabs(tz - patchz).gt.maxz).or.(dabs(tz-patchz).lt.minz)) then
	    i = i + s_natoms(targetsp)
	    cycle
	  end if
	  write(34,*) tz - patchz
	  ! Get bin position on surface
	  xbin = tx / xdelta
	  ybin = ty / ydelta
	  ! Store data
	  map(xbin,ybin,1) = map(xbin,ybin,1) + pvec(1)
	  map(xbin,ybin,2) = map(xbin,ybin,2) + pvec(2)
	  map(xbin,ybin,3) = map(xbin,ybin,3) + pvec(3)
	  map(xbin,ybin,4) = map(xbin,ybin,4) + (tz - patchz)
	  map(xbin,ybin,5) = map(xbin,ybin,5) + tz
	  !write(66,*) tz,patchz,dabs(tz - patchz)
	  acc(xbin,ybin) = acc(xbin,ybin) + 1.0
	  ! Accumulate global average
	  globalsum = globalsum + (tz - patchz)
	  globalsumpos = globalsumpos + tz
	  globalacc = globalacc + 1.0

	  i = i + s_natoms(targetsp)
	end do
	
	! Next frame (or finish)
116	if (nframes.EQ.endf) goto 801 
	goto 102

119	write(0,*) "HISTORY file ended prematurely!"
	write(0,"(A,I5,A)") "Managed ",nframesused," frames before error."

	goto 801

700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999

798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "Problem with HISTORY file - failed to read header."
	goto 999
800	write(0,*) "End of unformatted HISTORY file found."
801	write(0,*) "Finished."

	! Write out data
	open(unit=11,form="formatted",file=basename(1:baselen)//"map"//CHAR(48+targetsp)//".vf",status="replace")
	open(unit=12,form="formatted",file=basename(1:baselen)//"zmap"//CHAR(48+targetsp)//".surf",status="replace")
	open(unit=13,form="formatted",file=basename(1:baselen)//"zmap"//CHAR(48+targetsp)//".vf",status="replace")
	open(unit=14,form="formatted",file=basename(1:baselen)//"map"//CHAR(48+targetsp)//".surf",status="replace")
	open(unit=15,form="formatted",file=basename(1:baselen)//"zpos"//CHAR(48+targetsp)//".surf",status="replace")
	write(12,*) 2*NGRID+1, 2*NGRID+1
	write(12,"(6f8.4)") xdelta,0.0,0.0,0.0,ydelta,0.0
	write(12,"(3f10.4)") -0.5*cell(1),-0.5*cell(5),0.0
	write(12,*) "xyz"
	write(14,*) 2*NGRID+1, 2*NGRID+1
	write(14,"(6f8.4)") xdelta,0.0,0.0,0.0,ydelta,0.0
	write(14,"(3f10.4)") -0.5*cell(1),-0.5*cell(5),0.0
	write(14,*) "xyz"
	write(15,*) 2*NGRID+1, 2*NGRID+1
	write(15,"(6f8.4)") xdelta,0.0,0.0,0.0,ydelta,0.0
	write(15,"(3f10.4)") -0.5*cell(1),-0.5*cell(5),0.0
	write(15,*) "xyz"
	do ybin=-NGRID,NGRID
	  do xbin=-NGRID,NGRID
	    if (acc(xbin,ybin).gt.0.5) then
	      map(xbin,ybin,1) = map(xbin,ybin,1) / acc(xbin,ybin)
	      map(xbin,ybin,2) = map(xbin,ybin,2) / acc(xbin,ybin)
	      map(xbin,ybin,3) = map(xbin,ybin,3) / acc(xbin,ybin)
	      map(xbin,ybin,4) = map(xbin,ybin,4) / acc(xbin,ybin)
	      map(xbin,ybin,5) = map(xbin,ybin,5) / acc(xbin,ybin)
	    else
	      ! Set 'dead' points to the average surface value
	!write(0,*) globalsum,globalacc
              map(xbin,ybin,4) = globalsum / globalacc
              map(xbin,ybin,5) = globalsumpos / globalacc
	    end if
	    ! Average pointing vectors as function over surface
	    write(11,"(6f10.5)") xbin*xdelta, ybin*ydelta, 0.0, (map(xbin,ybin,n),n=1,3)
	    ! Surface zmap relative to calculation centre (patchz)
	    write(12,"(2f10.5)") map(xbin,ybin,4),acc(xbin,ybin)
	    ! Average pointing vectors as function over surface, including height component relative to average height
	    write(13,"(6f10.5)") xbin*xdelta, ybin*ydelta, map(xbin,ybin,5), (map(xbin,ybin,n),n=1,3)
	    ! Map of averages of z-components of vectors over surface
	    write(14,"(2f10.5)") map(xbin,ybin,3)
	    ! Height map of centres relative to average height
	    write(15,"(2f10.5)") map(xbin,ybin,5) - (globalsumpos / globalacc)
	  end do
	end do
	close(11)
	close(12)
	close(13)
	close(14)

	write(0,*) "Finished!"
999	close(10)
	close(20)

	end program gammadist
