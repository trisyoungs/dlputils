!	** glucgeom **
!	Output1: distances between each glucose's COM and a given chloride.
!	Output2: for each chloride, the ten closest OH distances (plus the owner glucose and the angles).

	program glucgeom
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,resfile
	character*20 :: temparg
	integer :: n,nframes,success,nargs,m,i,mol,aoff,cloff
	integer :: framestodo, nsites, t(9)
	integer :: species2(3), sp2 = 3
	integer, allocatable :: sites(:,:), cl_ownergluc(:), cl_glucsite(:)
	integer :: iargc
	real*8 :: a(3),b(3),c(3),temp(3),v1(3),v2(3),v1mag,v2mag,roo,angle,dp,maxoo,minang,totaverage,minoo
	real*8, allocatable :: cl_hdist(:), cl_hoangle(:), cl_comdist(:)
	real*8 :: calcangle, ncl

	minoo = 2.0

	nargs = iargc()
	if (nargs.NE.3) stop "Usage : glucgeom <HISTORYfile> <OUTPUTfile> <nframes> "
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temparg); read(temparg,"(I5)") framestodo
	
	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
        ! Now, read in the history header so that we have cell()
        if (readheader().EQ.-1) goto 799

	! Number of OH groups in molecule
	nsites = 5

	allocate(sites(nsites,2))
	allocate(cl_comdist(s_nmols(1)))
	allocate(cl_hdist(s_nmols(1)*nsites))
	allocate(cl_ownergluc(s_nmols(1)*nsites))
	allocate(cl_hoangle(s_nmols(1)*nsites))
	allocate(cl_glucsite(s_nmols(1)*nsites))

	sites(:,:) = reshape( (/ 8,9,10,11,12,20,21,22,23,24 /) , (/ nsites,2 /) )

	write(0,*) "OH site specification (O, H):"
	do m=1,nsites
	  write(0,*) (sites(m,n),n=1,2)
	end do
	write(0,"(A,A)") "Name of chloride species (check) : ",s_name(sp2)

	open(unit=18,file="cl-gluccom.dat",form="formatted",status="replace")
	open(unit=19,file="cl-hogeom.dat",form="formatted",status="replace")

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
101	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
	if (mod(nframes,100).EQ.0) write(0,*) nframes

	! Output 1: Chloride - glucose COM 
	write(18,"('FRAME ',I)") nframes
	call calc_com()
	i = s_start(sp2)-1
	do n=1,s_nmols(sp2)
	  c(1)=xpos(i+n)		! Chloride
	  c(2)=ypos(i+n)
	  c(3)=zpos(i+n)
	  do mol=1,s_nmols(1)
	    call pbc(comx(1,mol),comy(1,mol),comz(1,mol),c(1),c(2),c(3),temp(1),temp(2),temp(3))
	    cl_comdist(mol) = sqrt( (temp(1)-c(1))**2 + (temp(2)-c(2))**2 + (temp(3)-c(3))**2 )
	  end do
	  write(18,"(i3,1x,20f7.3)") n,cl_comdist
	end do

	! Output 2: 10 closest Cl-HO interactions (plus associated angles and owner glucose
	write(19,"('FRAME ',I)") nframes
	cloff = s_start(sp2)-1	! Atomic offset of chloride ions
	do n=1,s_nmols(sp2)		! Loop over chlorides

	  c(1)=xpos(cloff+n)		! Chloride position
	  c(2)=ypos(cloff+n)
	  c(3)=zpos(cloff+n)

	  cl_hdist = 0.0
	  cl_hoangle = 0.0
	  cl_ownergluc = 0
	  cl_glucsite = 0
	  i = 0

	  ! Fill cl_* arrays with distances, angles and parent glucose IDs

	  do mol=1,s_nmols(1)		! Loop over glucose molecules

	    aoff = (mol-1)*24		! Atomic offset of current glucose molecule

	    do m=1,nsites
	      i = i + 1
	      a(1)=xpos(sites(m,1)+aoff)	! Glucose Oxygen
	      a(2)=ypos(sites(m,1)+aoff)
	      a(3)=zpos(sites(m,1)+aoff)
	      b(1)=xpos(sites(m,2)+aoff)	! Glucose Hydrogen
	      b(2)=ypos(sites(m,2)+aoff)
	      b(3)=zpos(sites(m,2)+aoff)

	      ! gO-wO distance
	      call pbc(a(1),a(2),a(3),c(1),c(2),c(3),temp(1),temp(2),temp(3))
	      roo=sqrt( (temp(1)-c(1))**2 + (temp(2)-c(2))**2 + (temp(3)-c(3))**2 )
	      ! gO-gH-wO angle
              angle = calcangle(a,b,c)
	if (roo.lt.2.0) write(0,*) "Warning!!! Short Chloride-OH contact.",sites(m,1),cloff+n,nframes
	      cl_hdist(i) = roo
	      cl_hoangle(i) = angle
	      cl_ownergluc(i) = mol
	      cl_glucsite(i) = m

	    end do	! Loop over glucose OH sites

	  end do	! Loop over glucose molecules

	  !write(*,"(i3,1x,10(i2,1x,f7.3,f7.2,1x))") n,(cl_ownergluc(m),cl_hdist(m),cl_hoangle(m),m=80,100)

	  ! Sort the three arrays using the distance as the primary criteria
	  call sort4(s_nmols(1)*5,cl_hdist,cl_hoangle,cl_ownergluc,cl_glucsite)

	  ! Ouput the results (special case for isolated glucose)
	  if (s_nmols(1).eq.1) then
	    write(19,"(i3,1x,5(i2,1x,i1,1x,f7.3,f7.2,1x))") n,(cl_ownergluc(m),cl_glucsite(m),cl_hdist(m),cl_hoangle(m),m=1,5)
	  else
	    write(19,"(i3,1x,10(i2,1x,i1,1x,f7.3,f7.2,1x))") n,(cl_ownergluc(m),cl_glucsite(m),cl_hdist(m),cl_hoangle(m),m=1,10)
	  end if

	end do	! Loop over chlorides

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

	close(9)

	write(0,*) "Finished."
999	close(14)
	end program glucgeom

	real*8 function calcangle(a,b,c)
	use utility; implicit none
	real*8,dimension(3) :: a,b,c,v1,v2,temp1,temp2
	real*8 :: v1mag,v2mag,dp
        ! Angle between atoms a--b--c
	call pbc(a(1),a(2),a(3),b(1),b(2),b(3),temp1(1),temp1(2),temp1(3))
	call pbc(c(1),c(2),c(3),b(1),b(2),b(3),temp2(1),temp2(2),temp2(3))
	! Get the bond vectors
	v1 = b - temp1
	v2 = b - temp2
	v1mag = sqrt( v1(1)**2 + v1(2)**2 + v1(3)**2 )
	v2mag = sqrt( v2(1)**2 + v2(2)**2 + v2(3)**2 )
	v1 = v1 / v1mag
	v2 = v2 / v2mag
        ! Calculate dot product and angle...
        dp = (v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3))
        calcangle = acos(dp)* 57.29577951d0
	end function calcangle


      SUBROUTINE sort4(n,arr,brr,crr,drr)
      INTEGER n,M,NSTACK
      DOUBLE PRECISION arr(n),brr(n)
      INTEGER crr(n),drr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK),c,d
      DOUBLE PRECISION a,b,temp,tempi
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          c=crr(j)
          d=drr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
            crr(i+1)=crr(i)
            drr(i+1)=drr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
          crr(i+1)=c
          drr(i+1)=d
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        tempi=crr(k)
        crr(k)=crr(l+1)
        crr(l+1)=tempi
        tempi=drr(k)
        drr(k)=drr(l+1)
        drr(l+1)=tempi
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
          tempi=crr(l+1)
          crr(l+1)=crr(ir)
          crr(ir)=tempi
          tempi=drr(l+1)
          drr(l+1)=drr(ir)
          drr(ir)=tempi
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
          tempi=crr(l)
          crr(l)=crr(ir)
          crr(ir)=tempi
          tempi=drr(l)
          drr(l)=drr(ir)
          drr(ir)=tempi
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          temp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=temp
          tempi=crr(l+1)
          crr(l+1)=crr(l)
          crr(l)=tempi
          tempi=drr(l+1)
          drr(l+1)=drr(l)
          drr(l)=tempi
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
        c=crr(l)
        d=drr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        tempi=crr(i)
        crr(i)=crr(j)
        crr(j)=tempi
        tempi=drr(i)
        drr(i)=drr(j)
        drr(j)=tempi
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        crr(l)=crr(j)
        crr(j)=c
        drr(l)=drr(j)
        drr(j)=d
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
