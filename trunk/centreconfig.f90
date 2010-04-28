	! Takes the coordinates of a DL_POLY CONFIG file and outputs the
	! 'positive' version of them...
	
	! Variables
	character*100 discard,infile,outfile,header
	integer discardi
	real discardr
	real cell(9)
	character*5 :: atomname(400)
	integer image,info,atomno,n,genno,genflag,m,count
	integer :: iargc
	real*8 :: x(400),y(400),z(400),minx,miny,minz
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
	minx = 1e9; miny = 1e9; minz = 1e9

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
	! Check the image convention - if it is NOT zero we remove 3 lines
	IF (image.GT.0) THEN
	  READ(11,"(3F20.14)") cell
	END IF

	! Get the name of the output file....
	WRITE(*,*) "Name of the CONFIG file to output?"
	READ(*,*) outfile

	OPEN(UNIT=12,file=outfile,FORM="FORMATTED")
	WRITE(12,14) header
	WRITE(12,"(2I10)") 0,image
	IF (image.GT.0) THEN
	  WRITE(12,"(3F20.14)") cell
	END IF

	! Now do the atom position conversion
	count = 0
100	READ(11,17,END=200) atomname(count+1),discard,atomno
	IF (atomno.EQ.0) THEN
	  ! If we get here then the CONFIG is very young, i.e. no atom no's
	  ! So, we need to introduce a generic numbering system
	  genflag=1  ! Set the flag
	END IF
	! Read in the atomic coordinates
	count = count + 1
	READ(11,19) x(count),y(count),z(count)
	if (x(count).LT.0.0) x(count) = x(count) + cell(1)
	if (y(count).LT.0.0) y(count) = y(count) + cell(1)
	if (z(count).LT.0.0) z(count) = z(count) + cell(1)
	! Skip any other lines (forces/velocities)
	DO m=1,info
	  READ(11,14) discard
	END DO
	IF (genflag.EQ.1) THEN
	  atomno=n
	END IF
	goto 100

	! Adjust the coords as necessary....
200	do n=1,count
	  ! Now write out the next atomic data to the new CONFIG file....
	  WRITE(12,15) atomname(n),x(n),y(n),z(n)
	end do

999	CLOSE(11)
	close(12)
	END
