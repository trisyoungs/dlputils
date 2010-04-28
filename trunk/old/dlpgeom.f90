	! #########################################################
	! dlpgeom - analyses the average geometry in a history file
	! #########################################################

	program dlpgeom
	use utility; use dlprw
	implicit none
	real*8, dimension(:,:), allocatable :: gadata, gbdata, gtdata
	real*8, allocatable :: calctors(:)
	real*8  :: tx,ty,tz,vecji(3),vecjk(3),veckj(3),veckl(3),xp1(3),xp2(3),mag1,mag2,ktx,kty,ktz
        real*8  :: dp,angle,radcon,minimum,maximum,avg,sd,total
	character*80 :: geomfile,hisfile,temp
	character*6 :: geom
	integer :: nargs,n,m,i,j,k,l,status, nframes, framestodo, sp, a1, ngeom, ncount
	integer :: nbonds, nangles, ntorsions, success, nsp
	integer, allocatable :: bonds(:,:),angles(:,:),torsions(:,:),nspatoms(:),nmols(:),spstart(:)
	integer :: iargc
	logical :: finished = .false., calcbonds = .false., calcangles = .false., calctorsions = .false.
	real*8 :: dist, magnitude, dotproduct
1       FORMAT (A8/,3f20.14)
2       FORMAT (3F20.14)
3	FORMAT (2I10)

	radcon = 57.29577951d0

	nargs = iargc()
	if (nargs.NE.3) stop "Usage : dlpgeom <geom file> <history file> <nframes>"
	call getarg(1,geomfile)
	call getarg(2,hisfile)
	call getarg(3,temp)
	read(temp,"(I10)") framestodo

	! Read in the geometry specifications
	open(unit=10,file=geomfile,form='formatted',status="old")
	read(10,"(I3)") nsp
	write(0,*) nsp
	allocate(nspatoms(nsp)); allocate(nmols(nsp)); allocate(spstart(nsp))
	do n=1,nsp
	  read(10,"(I5)") nspatoms(n)
	write(0,*) nspatoms(n)
	  read(10,"(I5)") nmols(n)
	write(0,*) nmols(n)
	end do
	finished = .false.
5	read(10,"(A6,I10)") geom, ngeom
	select case (geom)
	  case ('bonds')
	    write(6,*) "Reading in list of bonds..."
	    calcbonds = .true.
	    nbonds = ngeom
	    allocate(bonds(nbonds,0:2))
	    do n=1,nbonds
	      read(10,"(3I4)") bonds(n,0),bonds(n,1),bonds(n,2)
	    end do
	    write(6,*) "Read in ",nbonds," bond definitions."
	  case ('angles')
	    write(6,*) "Reading in list of angles..."
	    calcangles = .true.
	    nangles = ngeom
	    allocate(angles(nangles,0:3))
	    do n=1,nangles
	      read(10,"(4I4)") angles(n,0),angles(n,1),angles(n,2),angles(n,3)
	    end do
	    write(6,*) "Read in ",nangles,"angle definitions."
	  case ('tors')
	    write(6,*) "Reading in list of torsions..."
	    calctorsions = .true.
	    ntorsions = ngeom
	    allocate(torsions(ntorsions,0:4))
	    allocate(calctors(ntorsions))
	    do n=1,ntorsions
	      read(10,"(5I4)") torsions(n,0),torsions(n,1),torsions(n,2),torsions(n,3),torsions(n,4)
	    end do
	    write(6,*) "Read in ",ntorsions,"torsion definitions."
	  case ('end')
	    finished = .true.
	  case default
	    stop "Geometry type not recognised!"
	end select
	if (.NOT.finished) goto 5
	close(10)

	! Calc natoms for comparison with history file
	i = 0
	do n=1,nsp
	  spstart(n) = i
	  i = i + nmols(n)*nspatoms(n)
	end do
	write(0,*) "Calculated",i," atoms from geometry file."

	! Open the history file
	call openhis(hisfile,15);
	success = readheader();
	if (success.NE.0) stop "Couldn't read history file header."
	if (natms.NE.i) stop "natoms in file differs from natoms specified in input!"

	allocate(gbdata(nbonds,framestodo))
	allocate(gadata(nangles,framestodo))
	allocate(gtdata(ntorsions,framestodo))
	gadata = 0.0
	gbdata = 0.0
	gtdata = 0.0

	nframes = 0
10	success = readframe()
	if (success.EQ.-1) then
	  write(0,*) "Error reading history file. Ending early."
	  goto 19
	end if
	if (success.EQ.1) then
	  write(0,*) "End of HISTORY file found."
	  goto 19
	end if
	nframes = nframes + 1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	if (calcbonds) then  ! Calculate distances...
	  do sp=1,nsp
	    a1 = spstart(sp)
	    do m=1,nmols(sp)
	      do n=1,nbonds
	        if (bonds(n,0).EQ.sp) then
	  	  i = bonds(n,1) + a1
	          j = bonds(n,2) + a1
		  call pbc(xpos(i),ypos(i),zpos(i),xpos(j),ypos(j),zpos(j),tx,ty,tz)
	          gbdata(n,nframes) = gbdata(n,nframes) + dist(tx,ty,tz,j)
	      ! Write out distance data HACK
	      ! write(88,"(i6,i6,10f12.4)") nframes,m,calctors(n)
	      write(86,"(10f12.4)") dist(tx,ty,tz,j)
		end if
	      end do
	      a1 = a1 + nspatoms(sp)
	    end do
	  end do
	end if

	if (calcangles) then  ! Calculate angles....
	  do sp=1,nsp
	    a1 = spstart(sp)
	    do m=1,nmols(sp)
	      do n=1,nangles
		if (angles(n,0).EQ.sp) then
		  ! Minimum Image w.r.t. atom j ('elbow')
		  j=angles(n,2) + a1
		  ! Vector j->i
		  i=angles(n,1) + a1
		  call pbc(xpos(i),ypos(i),zpos(i),xpos(j),ypos(j),zpos(j),tx,ty,tz)
		  call getvector(tx,ty,tz,xpos(j),ypos(j),zpos(j),vecji(1),vecji(3),vecji(3))
		  mag1 = magnitude(vecji)
		  ! Vector j->k
		  k=angles(n,3) + a1
		  call pbc(xpos(k),ypos(k),zpos(k),xpos(j),ypos(j),zpos(j),tx,ty,tz)
		  call getvector(tx,ty,tz,xpos(j),ypos(j),zpos(j),vecjk(1),vecjk(3),vecjk(3))
		  mag2 = magnitude(vecjk)

		  ! Calculate dot product and angle...
		  dp = dotproduct(vecji, vecjk) / (mag1 * mag2)
		  angle = acos(dp)*radcon
		  gadata(n,nframes) = gadata(n,nframes) + angle
		end if
	      end do
	      a1 = a1 + nspatoms(sp)
	    end do
          end do
	end if

	if (calctorsions) then  ! Calculate torsions....
	  do sp=1,nsp
	    a1 = spstart(sp)
	    do m=1,nmols(sp)
	      do n=1,ntorsions
		if (torsions(n,0).EQ.sp) then
		  ! Minimum Image w.r.t. atom j 
		  j=torsions(n,2) + a1
		  ! Angle i-j-k
		  ! Vector j->i
		  i=torsions(n,1) + a1
		  call pbc(xpos(i),ypos(i),zpos(i),xpos(j),ypos(j),zpos(j),tx,ty,tz)
		  call getvector(tx,ty,tz,xpos(j),ypos(j),zpos(j),vecji(1),vecji(2),vecji(3))
		  ! Vector j->k
		  k=torsions(n,3) + a1
		  call pbc(xpos(k),ypos(k),zpos(k),xpos(j),ypos(j),zpos(j),ktx,kty,ktz)
		  call getvector(ktx,kty,ktz,xpos(j),ypos(j),zpos(j),vecjk(1),vecjk(2),vecjk(3))
		  ! Angle j-k-l (mim w.r.t. k (mim j))
		  ! Vector k->j
		  veckj = -vecjk
		  ! Vector k->l
		  l=torsions(n,4) + a1
		  call pbc(xpos(l),ypos(l),zpos(l),ktx,kty,ktz,tx,ty,tz)
		  call getvector(tx,ty,tz,ktx,kty,ktz,veckl(1),veckl(2),veckl(3))

		  ! Calculate cross products and magnitudes
		  call crossproduct(vecji,vecjk,xp1)
		  mag1 = magnitude(xp1)
		  call crossproduct(veckj,veckl,xp2)
		  mag2 = magnitude(xp2)
		  
		  ! Calculate dot product and angle...
		  dp = dotproduct(xp1, xp2) / (mag1 * mag2)
		  angle = acos(dp)*radcon

		  ! Calculate sign
		  dp = dotproduct(xp1, veckl)
		  if (dp.lt.0) angle = -angle

		  calctors(n) = angle
	      ! Write out torsion data HACK
	      ! write(88,"(i6,i6,10f12.4)") nframes,m,calctors(n)
	      write(88,"(10f12.4)") calctors(n)
		  gtdata(n,nframes) = gtdata(n,nframes) + angle
		end if
	      end do
	      a1 = a1 + nspatoms(sp)
	    end do
          end do
	  ! write(88,"(I10,1x,20(F5.1,1x))") nframes,(gtdata(n,nframes),n=1,ntorsions)
	end if

	if (nframes.EQ.framestodo) goto 19
	goto 10


19	write(0,"(A)") "Analysing data...."
	close(14)
20	FORMAT ('Bond    ',I2,' : Atoms ',6x,2I3,4(2x,f12.6))
21	FORMAT ('Angle   ',I2,' : Atoms ',3x,3I3,4(2x,f12.6))
22	FORMAT ('Torsion ',I2,' : Atoms ',4I3,4(2x,f12.6))
23	FORMAT ('                                   Min           Max           Avg          S.D.')

	! Normalise w.r.t. number of molecules
	do n=1,nframes
	  do m=1,nbonds
	    gbdata(m,n) = gbdata(m,n) / nmols(bonds(m,0))
	  end do
	  do m=1,nangles
	    gadata(m,n) = gadata(m,n) / nmols(angles(m,0))
	  end do
	  do m=1,ntorsions
	    gtdata(m,n) = gtdata(m,n) / nmols(torsions(m,0))
	  end do
	end do

	! Calculate statistics
	write(6,23)
	do n=1,nbonds
	  minimum = minval(gbdata(n,:))
	  maximum = maxval(gbdata(n,:))
	  avg = sum(gbdata(n,:)) / real(nframes)
	  ! Now for S.D.
	  total = 0.0d0
	  do m=1,nframes
	    total = total + (gbdata(n,m) - avg)**2
	  end do
	  sd = SQRT( total / real(nframes) )
	  write(6,20) n,(bonds(n,m),m=1,2),minimum,maximum,avg,sd
	end do

	do n=1,nangles
	  minimum = minval(gadata(n,:))
	  maximum = maxval(gadata(n,:))
	  avg = sum(gadata(n,:)) / real(nframes)
          ! Now for S.D.
          total = 0.0d0
          do m=1,nframes
            total = total + (gadata(n,m) - avg)**2
          end do
          sd = SQRT( total / real(nframes) )
          write(6,21) n,(angles(n,m),m=1,3),minimum,maximum,avg,sd
        end do

	do n=1,ntorsions
	  minimum = 200.0; maximum = -200.0; avg = 0.0
	  do m=1,nframes
	    if (gtdata(n,m).LT.minimum) minimum = gtdata(n,m)
	    if (gtdata(n,m).GT.maximum) maximum = gtdata(n,m)
	    avg = avg + gtdata(n,m)
	  end do
	  avg = avg / real(nframes)
	  !minimum = minval(gtdata(n,:))
	  !maximum = maxval(gtdata(n,:))
	  !avg = sum(gtdata(n,:)) / real(nframes)
          ! Now for S.D.
          total = 0.0d0
          do m=1,nframes
            total = total + (gtdata(n,m) - avg)**2
          end do
          sd = SQRT( total / real(nframes) )
          write(6,22) n,(torsions(n,m),m=1,4),minimum,maximum,avg,sd
        end do

	end program dlpgeom
	
        real*8 function dist(ix,iy,iz,j)
        use dlprw; implicit none
        integer :: j
        real*8 :: ix,iy,iz
        ! Calculate the distance between the two atoms a (specified by ix,iy,iz) and j.
        dist = SQRT( (ix-xpos(j))**2 + (iy-ypos(j))**2 + (iz-zpos(j))**2 )
        end function dist

        real*8 function magnitude(vec)
	implicit none
        real*8 :: vec(3)
        magnitude=dsqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
        end function magnitude

        subroutine crossproduct(abc,xyz,result)
	implicit none
        real*8 :: abc(3), xyz(3), result(3)
        result(1)=abc(2)*xyz(3) - abc(3)*xyz(2)
        result(2)=abc(3)*xyz(1) - abc(1)*xyz(3)
        result(3)=abc(1)*xyz(2) - abc(2)*xyz(1)
        end subroutine crossproduct

        real*8 function dotproduct(abc,xyz)
	implicit none
        real*8 :: abc(3), xyz(3), result
        result=abc(1)*xyz(1) + abc(2)*xyz(2) + abc(3)*xyz(3)
        dotproduct = result
        end function dotproduct

