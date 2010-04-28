!	Converts a BGF file of a system exported from Cerius2.

	IMPLICIT NONE
!	Variables
	CHARACTER*100 dis
	CHARACTER*80 simname,bgffile,configfile
	CHARACTER*5 atomname(8000),temp(4),vdw(50,2)
	CHARACTER*5 uniqueatoms(50)
	INTEGER disn,n,m,o,p,numatoms,numcon
	INTEGER con(12)
	INTEGER pcond,crystal
	INTEGER species(10,2)       ! (x,1)=natoms in species x, (x,2)=nummols of x
	REAL x(8000),y(8000),z(8000),disr,cell(3),charge(8000),mass(8000)
	INTEGER found,terminate,ligfreeze,specifiedatoms

!	Formats
10	FORMAT (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,f10.5,f10.5,f10.5,1x,a5,i3,i2,1x,f8.5)
11      FORMAT (A6,12I6)
12	FORMAT (A5,4X,A5,2F12.2)
13	FORMAT (3A5,2F10.2)
14	FORMAT (4A5,2F12.2,I7)
15	FORMAT (A5/,3f20.14)
16	FORMAT (3F20.14)
17	FORMAT (A5,3X,A5,4X,A2,1X,2F12.5)

!	Preprocess the bgf file
	WRITE (*,*) "Name of bgf file to convert:"
	READ (*,*) bgffile
	OPEN (ERR=991,UNIT=12,FILE=bgffile,FORM="FORMATTED")
!       Get the number of atoms and connection lines in the file
50	READ(12,"(A100)") dis
	IF (dis(:6).NE."FORMAT") GOTO 50
	numatoms=0
60      READ(12,"(A100)") dis
        numatoms=numatoms+1
        IF (dis(:6).NE."FORMAT") GOTO 60
        numatoms=numatoms-1
	WRITE(*,*) "Atoms in the system : ",numatoms
        REWIND(12)
	! See if we can find crystal dimensions....
	crystal=0
70	READ(12,"(A100)") dis
	IF (dis(:6).EQ."CRYSTX") THEN
	  ! Found crystal data - parse the cell dimensions....
	  WRITE(*,*) "Got crystal dimensions from bgf file..."
	  READ(dis(11:18),"(F8.5)") cell(1)
	  READ(dis(22:29),"(F8.5)") cell(2)
	  READ(dis(33:40),"(F8.5)") cell(3)
	  crystal=1
	END IF
	IF (dis(:6).NE."FORMAT") GOTO 70

!	Read in the list of atoms:names, coordinates e.t.c.
	REWIND(12)
80	READ(12,"(A100)") dis
        IF (dis(:6).NE."FORMAT") GOTO 80
	DO n=1,numatoms
	  READ(12,10) dis,disn,dis,dis,dis,dis,x(n),y(n),z(n),atomname(n),disn,disr,charge(n)
	END DO

!	Now we can write the CONFIG and FIELD files!!!!
	WRITE(*,*) "Name of CONFIG file to output: (CONFIG)"
	READ(*,*) configfile
	IF (configfile.EQ."") configfile="CONFIG"
	WRITE(*,*) "Name of the simulation?"
	READ(*,*) simname
	IF (simname.EQ."") simname="***Untitled"
!	Ask which periodic conditions to use...
	WRITE(*,*) "Periodic conditions - 0=none ,1=cubic, 2=orthorhombic"
	READ(*,"(I4)") pcond
	IF (crystal.EQ.1) THEN
	  WRITE(*,*) "Using box dimensions from bgf file..."
	ELSE
	  IF (pcond.EQ.1) THEN
	    WRITE(*,*) "-- Enter length of box side (MUST be a *REAL* number)"
	    READ(*,"(F10.4)") cell(1)
	    cell(2)=cell(1)
	    cell(3)=cell(1)
	  END IF
	  IF (pcond.EQ.2) THEN
	    WRITE(*,*) "-- Enter lengths of box sides (in REAL format) - Return after each one:"
	    READ(*,"(F10.4)") cell(1)
	    READ(*,"(F10.4)") cell(2)
	    READ(*,"(F10.4)") cell(3)
	  END IF
	END IF
        OPEN (ERR=992,UNIT=13,FILE=configfile, FORM="FORMATTED")
	WRITE(13,*) simname
	IF (pcond.GT.0) THEN
	  WRITE(13,"(2I10)") 0,pcond
	  WRITE(13,16) cell(1),0.0,0.0,0.0,cell(2),0.0,0.0,0.0,cell(3)
	ELSE
	  WRITE(13,"(2I10)") 0,0
	END IF
	DO n=1,numatoms
	  WRITE(13,15) atomname(n),x(n),y(n),z(n)
	END DO
	CLOSE(13)
	GOTO 998

991	WRITE(*,*) "ERROR: BGF file not found."
	STOP
992	WRITE(*,*) "ERROR: Couldn't open CONFIG file."
	STOP
998	WRITE(*,*) "Finished."
999	END
