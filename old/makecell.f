!	Creates an NxN replicated FCC cell.

	IMPLICIT NONE
!	Variables
	CHARACTER*100 dis
	CHARACTER*80 simname,bgffile,configfile
	CHARACTER*5 atomname(2)
	CHARACTER*5 uniqueatoms(50)
	integer n,m,i,j,k,keytrj,imcon
	real*8 atomx(8),atomy(8),atomz(8),cell(9)

!	Formats
10	FORMAT (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,f10.5,f10.5,f10.5,1x,a5,i3,i2,1x,f8.5)
11      FORMAT (A6,12I6)
12	FORMAT (A5,4X,A5,2F12.2)
13	FORMAT (3A5,2F10.2)
14	FORMAT (4A5,2F12.2,I7)
15	FORMAT (A5/,3f20.14)
16	FORMAT (3F20.14)
17	FORMAT (A5,3X,A5,4X,A2,1X,2F12.5)

	data atomx/0.0d0,0.5d0,0.0d0,0.5d0,0.0d0,0.0d0,0.5d0,0.5d0/
	data atomy/0.0d0,0.5d0,0.5d0,0.0d0,0.0d0,0.5d0,0.0d0,0.5d0/
	data atomz/0.0d0,0.0d0,0.5d0,0.5d0,0.5d0,0.0d0,0.0d0,0.5d0/
	data cell/1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0/

	atomname(1)="P    "
	atomname(2)="N    "

	WRITE(*,*) "Replications in each direction?"
	READ(*,*) M

!	Now we can write the CONFIG and FIELD files!!!!
	WRITE(*,*) "Name of CONFIG file to output: (CONFIG)"
	READ(*,*) configfile
	IF (configfile.EQ."") configfile="CONFIG"
	simname="FCC lattice ("//CHAR(48+M)//"x"//CHAR(48+M)//")"
	imcon=1
	do n=1,9
	  cell(n)=cell(n)*real(M)
	end do
	keytrj=0
	
        OPEN (ERR=992,UNIT=13,FILE=configfile, FORM="FORMATTED")
	WRITE(13,*) simname
	WRITE(13,"(2I10)") keytrj,imcon
	WRITE(13,16) cell
	do i=0,M-1
	  do j=0,M-1
	    do k=0,M-1
	      do n=1,8
		if (n.LE.4) write(13,"(A5)") atomname(1)
	        if (n.GE.5) write(13,"(A5)") atomname(2)
		write(13,16) atomx(n)+real(i),atomy(n)+real(j),atomz(n)+real(k)
	      end do
	    end do
	  end do
	end do

	CLOSE(13)
	GOTO 998

992	WRITE(*,*) "ERROR: Couldn't open CONFIG file."
	STOP
998	WRITE(*,*) "Finished."
999	END
