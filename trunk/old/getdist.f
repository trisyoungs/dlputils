	! Take a config file and get the given distance (over all occurrences).
	EXTERNAL outinfo
	! Variables
	integer maxs
	PARAMETER(MAXS=20)
	character*100 discard,infile,outfile,header,dlpoutfile
	integer discardi,outinfo
	real discardr
	real*8 cell(9),x,y,z,fx,fy,fz,vx,vy,vz
	real*8 xx(50),yy(50),zz(50)
	real*8 minval,maxval,average
	character*5 atomname(10,50),t1,t2
	integer image,info,atomno,n,genno,genflag,m
	integer edittype,editatom,o,p,atom1,atom2,count
	character*20 s_name(MAXS)
	integer nspecies,s_natoms(MAXS),s_nmols(MAXS),s_start(MAXS)
	COMMON /species/ s_name,s_natoms,s_nmols,s_start,nspecies

	! Formats
	! 'Discard' line / DL_POLY header line
14	FORMAT (A100)
15	FORMAT (A5/,3f20.14)
	! DL_Poly CONFIG 1 - data in file/image convention/???
16	FORMAT (I10,I10,A100)
	! DL_POLY CONFIG 2 - atom name/xxx/no
17	FORMAT (A5,A12,I3)
	! DL_POLY CONFIG 3 - atom positions OR forces OR velocities
18	FORMAT (A3,F17.14,A3,F17.14,A3,F17.14)
19	FORMAT (3F20.14)

	! Init
	genflag=0

	! Get the file name of the CONFIG to work on....
	WRITE(*,*) "Name of the DL_POLY CONFIG file?"
	READ(*,*) infile
	
	OPEN(UNIT=11,file=infile,FORM="FORMATTED")

	! First determine which data are in the CONFIG file
	WRITE(*,*) "CONFIG file data:"
	READ(11,14) header      ! Remove the first line (simulation name)
	READ(11,16) info,image,discard
	IF (info.EQ.0) WRITE(*,*) "Type : Coordinates only."
	IF (info.EQ.1) WRITE(*,*) "Type : Coordinates and velocities."
	IF (info.EQ.2) WRITE(*,*) "Type : Coordinates, velocities and forces."
	! Check the image convention - if it is NOT zero we read the cell dimensions
	IF (image.GT.0) THEN
	  READ(11,"(3F20.14)") cell
	END IF

	! Get the molecule data from the user....
50	WRITE(0,*) "Name of the outfile to pull molecule data from? (0 for none)"
	READ(*,*) dlpoutfile
	IF (dlpoutfile(1:2).EQ."0 ") THEN  ! Get manual input.
	  WRITE(*,*) "Number of molecule types in config file?"
	  READ(*,*) nspecies
	  DO n=1,nspecies
	    WRITE(*,*) "Number of molecules of type: ",n,"?"
	    READ(*,*) s_nmols(n)
	    WRITE(*,*) "Number of atoms in molecule type: ",n,"?"
	    READ(*,*) s_natoms(n)
	  END DO
	ELSE     ! Read the data from an OUT file....
	  OPEN(UNIT=13,FILE=dlpoutfile,FORM="FORMATTED",IOSTAT=n)
	  IF (n.EQ.0) THEN   ! File opened successfully....
	    IF (outinfo(13,1).EQ.-1) GOTO 988
	    success=0
	  ELSE
	    success=-1
	  END IF
	END IF
	IF (success.EQ.-1) GOTO 50

	! Now, read in the 'old' atom names from the CONFIG file....
	WRITE(*,*) "Processing CONFIG file...."
	DO n=1,nspecies
	 DO o=1,s_nmols(n)
	  DO m=1,s_natoms(n)
	    READ(11,17,END=200,ERR=200) atomname(n,m)
	    ! Skip other lines.....
	    DO p=0,info
	      READ(11,14,END=200,ERR=200) discard
	    END DO
	  END DO
	 END DO
	END DO

	! Print the list of atom names now and ask for the target....
70	WRITE(*,*) ""
	DO n=1,nspecies
	  WRITE(*,*) "Species ",n,":"
	  WRITE(*,"(8(I3,1X,'(',A5,')',3x))") (INT(X),atomname(n,X),X=1,s_natoms(n))
	END DO
	WRITE(*,*) ""
	WRITE(*,*) "Select species containing the atoms:"
	READ(*,*) edittype
	WRITE(*,*) " Enter the IDs of the two atoms:"
	READ(*,*) atom1
	READ(*,*) atom2

	! Rewind and prepare the CONFIG file...
90	REWIND(11)
	READ(11,14) header
	READ(11,16) info,image,discard
	IF (image.GT.0) THEN
	  READ(11,"(3F20.14)") cell
	END IF

	minval=9999.9
	maxval=0.0
	average=0.0
	count=0
	! 'Locate' the correct species in the CONFIG file.....
	DO n=1,nspecies
	  IF (edittype.EQ.n) THEN    ! This is the target species...
	    DO m=1,s_nmols(n)
	      DO o=1,s_natoms(n)
	        READ(11,17) t1,discard,atomno   ! Don't really care about this...	
	        READ(11,19) xx(o),yy(o),zz(o)
	        IF (info.GT.0) READ(11,19) vx,vy,vz
	        IF (info.GT.1) READ(11,19) fx,fy,fz
	      END DO
	      ! Now have the coordinates of the whole species, so calc distance and store..
	      ! Perform minimum image checks....
	      If (xx(atom1)-xx(atom2).GT.cell(1)/2.0) xx(atom2)=xx(atom2)+cell(1)
	      If (xx(atom1)-xx(atom2).LT.-cell(1)/2.0) xx(atom2)=xx(atom2)-cell(1)
	      If (yy(atom1)-yy(atom2).GT.cell(5)/2.0) yy(atom2)=yy(atom2)+cell(5)
	      If (yy(atom1)-yy(atom2).LT.-cell(5)/2.0) yy(atom2)=yy(atom2)-cell(5)
	      If (zz(atom1)-zz(atom2).GT.cell(9)/2.0) zz(atom2)=zz(atom2)+cell(9)
	      If (zz(atom1)-zz(atom2).LT.-cell(9)/2.0) zz(atom2)=zz(atom2)-cell(9)
	      dx=xx(atom1)-xx(atom2)
	      dy=yy(atom1)-yy(atom2)
	      dz=zz(atom1)-zz(atom2)
	      dist=SQRT(dx**2 + dy**2 + dz**2)
	      IF (dist.LT.minval) minval=dist
	      IF (dist.GT.maxval) maxval=dist
	      average=average+dist
	      count=count+1
	    END DO
	  ELSE
	    ! Just read in the atoms and continue...
	    DO m=1,s_nmols(n)
	      DO o=1,s_natoms(n)
	        READ(11,17) t1,discard,atomno
	        READ(11,19) x,y,z
	        IF (info.GT.0) READ(11,19) vx,vy,vz
	        IF (info.GT.1) READ(11,19) fx,fy,fz
	      END DO
	    END DO
	  END IF
	END DO

	! Write out the results...
	WRITE(*,"('Average distance between atom pair is ',F8.6,' Angstroms (',I4,' occurrences).')") average/count,count
	WRITE(*,"('Minimum value found was ',F8.6,' and maximum was ',F8.6,' Angstroms.')") minval, maxval
	WRITE(*,*) "Enter another atom pair or '0' for a different species:   (-1 to exit)"
	READ(*,*)  atom1
	IF (atom1.EQ.0) GOTO 70
	IF (atom1.EQ.-1) GOTO 990
	READ(*,*) atom2
	GOTO 90

	! Exit messages
200	WRITE(*,*) "Error while reading in CONFIG file data."
	WRITE(*,*) "Were the molecular species specified correctly?"
	DO p=1,nspecies
	  WRITE(*,*) "  Type ",p,": nmols=",s_nmols(p),"  natoms=",s_natoms(p)
	END DO
	WRITE(*,*) "Error occured on Type(",n,") Mol(",o,") Atom(",m,")."
	GOTO 999

988	WRITE(0,*) "Error while reading OUT file."
	GOTO 999
990	WRITE(*,*) "Finished!"
999	CLOSE(11)
	CLOSE(12)
	END
