!	** addcomatom **
!	Adds in an extra atom to the config file for a species, with coords = COM
	IMPLICIT NONE
	integer maxn,maxframes,maxbins,readframe,maxs,outinfo,gcatomlist(50),patoms(2)
	real*8 mimdist
	PARAMETER(MAXN=8000,MAXFRAMES=1100,MAXBINS=400,MAXS=20)
	character*80 header,infile,outfile,dlpoutfile,discard80
	character*8 atmnam,speciesn(2),comname
	character*8 discard
	integer success,calctype,species(2),speciesa(2),bin,X,target,counter,comlimit
	integer keytrj, imcon, natms, nstep, iatm, n, discardn, nframes, l, fileform, m, o, p
	real*8 tstep, binwidth, dist
	real*8 cell(9), xpos(MAXN),ypos(MAXN),zpos(MAXN), charge(MAXN), mass(MAXN)
	real*8 tx,ty,tz,discardr,c1x,c1y,c1z,c2x,c2y,c2z
	real*8 s_mass(MAXS,20),comx,comy,comz
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS),centre(2),comtype

14      FORMAT (A80)
15      FORMAT (A5/,3f20.14)
	! DL_Poly CONFIG 1 - data in file/image convention/???
16      FORMAT (I10,I10,A80)
	! DL_POLY CONFIG 2 - atom name/xxx/no
17      FORMAT (A5,A12,I3)
	! DL_POLY CONFIG 3 - atom positions OR forces OR velocities
18      FORMAT (A3,F17.14,A3,F17.14,A3,F17.14)
19      FORMAT (3F20.14)

	! Type of centre-of-mass calculations...
	comtype=1	     ! 1 = centre of mass,     2 = geometric centre of mass
	WRITE(0,*) "*** rdf"
	!IF (comtype.EQ.1) WRITE(0,*) "( using true centres-of-mass for calculations )"
	!IF (comtype.EQ.2) WRITE(0,*) "( using geometric centres for calculations )"

	! Set the parameters for the calculation........
	nspecies=2
	s_nmols(1)=200
	s_nmols(2)=200
	s_natoms(1)=10
	s_natoms(2)=1
	comlimit=8
	s_mass(1,1)=14.0067
	s_mass(1,2)=12.0107
	s_mass(1,3)=14.0067
	s_mass(1,4)=12.0107
	s_mass(1,5)=12.0107
	s_mass(1,6)=1.0079
	s_mass(1,7)=1.0079
	s_mass(1,8)=1.0079
	s_mass(1,9)=15.0344
	s_mass(1,10)=15.0344
	s_mass(2,1)=35.453
	comname="Chg"
	target=1

	natms=0
	DO n=1,nspecies
	  natms=natms+s_nmols(n)*s_natoms(n)
	END DO
	WRITE(0,*) "(Expecting ",natms," atoms in CONFIG file)" 

	WRITE(*,*) "Name of CONFIG file?"
	READ(*,*) infile
	OPEN(UNIT=11,FILE=infile,FORM="FORMATTED")
	OPEN(UNIT=12,FILE="out.CONFIG",FORM="FORMATTED")

	! Now, must construct the position and mass arrays...
	READ(11,14) header      ! Remove the first line (simulation name)
	READ(11,16) keytrj,imcon,discard80
	IF (keytrj.EQ.0) WRITE(*,*) "Type : Coordinates only."
	IF (keytrj.EQ.1) WRITE(*,*) "Type : Coordinates and velocities."
	IF (keytrj.EQ.2) WRITE(*,*) "Type : Coordinates, velocities and forces."
	! Check the image convention - if it is NOT zero we remove 3 lines
	IF (imcon.GT.0) THEN
	  READ(11,"(3F20.14)") cell
	END IF

	! Write out the header for the new config file...
	WRITE(12,14) header
	WRITE(12,"(2I10)") 0,imcon
	IF (imcon.GT.0) THEN
	  WRITE(12,"(3F20.14)") cell
	END IF

	! Main loop begins here. Do all 'on the fly'

	DO n=1,nspecies
	DO m=1,s_nmols(n)
	
	! Read in the coords of this molecule
	DO o=1,s_natoms(n)
	  READ(11,17) atmnam,discard,discard80
	  READ(11,19) xpos(o),ypos(o),zpos(o)
	  ! Skip any other lines (forces/velocities)
	  DO p=1,keytrj
	    READ(11,14) discard
	  END DO
	  ! Can write them out again straight away....
	  WRITE(12,15) atmnam,xpos(o),ypos(o),zpos(o)
	END DO

	! Can write them out straight away....

	! If this is the target species, calc COM for new atom...
	IF (n.EQ.target) THEN
	  ! Calculate COM
	  ! Add the first atom....
	  comx=xpos(1)
	  comy=ypos(1)
	  comz=zpos(1)
	  DO p=2,comlimit   
	    CALL PBC(xpos(p),ypos(p),zpos(p),comx,comy,comz,tx,ty,tz,cell)
	    ! Add to the GCpos : multiply up first, then add, then divide for new GC
	    comx = ( (comx*(p-1)) + tx) / p
	    comy = ( (comy*(p-1)) + ty) / p
	    comz = ( (comz*(p-1)) + tz) / p
	  END DO
	  WRITE(12,15) comname,comx,comy,comz
	END IF
	
	END DO
	END DO

	GOTO 801

700	WRITE(*,*) "INFILE and OUTFILE have the same name!!!"
	GOTO 999

798	WRITE(0,*) "Problem with OUTPUT file."
	GOTO 999
799	WRITE(0,*) "HISTORY file ended prematurely!"
	GOTO 801
800	WRITE(0,*) "End of unformatted HISTORY file found."

	! Write out the new CONFIG file.......
801	WRITE(0,*) "Finished"

999	CLOSE(10)
	CLOSE(13)
	END

	SUBROUTINE PBC(x1,y1,z1,x2,y2,z2,x3,y3,z3,cell)
	! Performs minimum image convention.
	! xyz1: point to consider; xyz2: reference point; xyz3: minimum image coords of xyz1 (result)
	REAL*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,cell(9)
	! Performs minimum image convention.....
	x3=x1 - cell(1)*NINT((x1-x2)/cell(1))
	y3=y1 - cell(5)*NINT((y1-y2)/cell(5))
	z3=z1 - cell(9)*NINT((z1-z2)/cell(9))
	END


