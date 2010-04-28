!	Test program
	EXTERNAL sub
	INTEGER n,m


	v1x=1.0/SQRT(2.0)
	v1y=1.0/SQRT(2.0)
	v1z=0.0
	v2x=0.0
	v2y=1.0
	v2z=0.0
	
	WRITE(0,*) ACOS(v1x*v2x + v1y*v2y + v1z*v2z)
	WRITE(0,*) (ACOS(v1x*v2x + v1y*v2y + v1z*v2z) / 3.14159) * 180.0
	n=15
	m=6
	WRITE(0,*) n/m
	WRITE(0,*) MOD(n,m)
	n=5
	m=101
	WRITE(0,*) " n set to ",n

	CALL sub(n)

	WRITE(0,*) " n is now ",n," in main program"
	END

	SUBROUTINE sub(a)
	INTEGER a
	WRITE(0,*) m
	a=a+5
	WRITE(0,*) "Inside subroutine, n is ",a," after addition"
	END
