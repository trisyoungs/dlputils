	program ortho
	real*8 :: i(3), j(3), k(3), vx, vy, vz, temp(3)
	real*8 :: ao(3), ax(3), ay(3), az(3), dp, xmag, ymag, zmag

	!i(:) = (/ -1.16112000000000 , 1.18815000000000 , 0.34488000000000 /)
	!j(:) = (/ -1.92586000000000 , 2.23706000000000 , -0.1731400000000 /)
	!k(:) = (/ -2.83079000000000 , 1.88064000000000 , -0.2597700000000 /)
	i(:) = (/ -1.16112000000000 ,  1.18815000000000  , 0.34488000000000 /)
	j(:) = (/ 1.07981000000000 ,  0.53609000000000  , 0.59846000000000 /)
	k(:) = (/ 0.18479000000000 ,  1.48753000000000 ,  0.01872000000000 /)
	!i(:) = (/ 0.0, 0.0, 0.0 /)
	!j(:) = (/ 1.0,-1.0, 1.0 /)
	!k(:) = (/ 1.0, 0.0, 1.0 /)

	! X axis : i->j
	call getvector(j(1),j(2),j(3),i(1),i(2),i(3),vx,vy,vz)
	write(0,"(A,3F10.6)") "vecX  ",vx,vy,vz
	xmag=SQRT(vx**2 + vy**2 + vz**2)     ! Magnitude of vector
	ax(1)=vx 	! Store the un-normalised components of the x-axis vector
	ax(2)=vy 
	ax(3)=vz 

	! Set the origin of the axis system - halfway along the x-axis vector
	ao(1)=i(1) + 0.5*vx
	ao(2)=i(2) + 0.5*vy
	ao(3)=i(3) + 0.5*vz

	! Y-axis : ao->k
	call getvector(k(1),k(2),k(3),i(1),i(2),i(3),vx,vy,vz)
	write(0,"(A,3F10.6)") "vecY  ",vx,vy,vz
	ay(1)=vx	! Store the un-normalised components of the y-axis vector
	ay(2)=vy
	ay(3)=vz

	! Orthogonalise this vector w.r.t. ax()
	write(0,"(A,3F10.6)") "vecYN ",ay
	dp = ax(1)*ay(1) + ax(2)*ay(2) + ax(3)*ay(3)
	write(0,"(A,F10.6)") "DP    ",dp
	temp = ay - (dp / xmag**2)*ax
	ay = temp

	! Normalise the X and Y vectors
	ymag=SQRT(ay(1)**2 + ay(2)**2 + ay(3)**2)     ! Magnitude of vector
	write(0,"(A,3F10.6)") "vecY O",ay
	ay=ay/ymag  
	ax=ax/xmag
	write(0,"(A,3F10.6)") "vecYNO",ay
	write(0,"(A,3F10.6)") "vecXNO",ax

	! 4) Calculate the z axis from the cross product of the x and y axes....
	az(1)=ax(2)*ay(3) - ax(3)*ay(2)
	az(2)=ax(3)*ay(1) - ax(1)*ay(3)
	az(3)=ax(1)*ay(2) - ax(2)*ay(1)

	write(0,"(A3,3F20.14)") "X: ",ax(1),ax(2),ax(3)
	write(0,"(A3,3F20.14)") "Y: ",ay(1),ay(2),ay(3)
	write(0,"(A3,3F20.14)") "Z: ",az(1),az(2),az(3)

	open(unit=9,file="y.CONFIG",form="formatted",status="replace")
	write(9,"(A)") "Original MSI: betaglucosecell_4.msi"
	write(9,"(2I10)") 0,1
	write(9,"(3F20.14)") 50.0,0.0,0.0
	write(9,"(3F20.14)") 0.0,50.0,0.0
	write(9,"(3F20.14)") 0.0,0.0,50.0
	write(9,"('N',/,3F20.14)") 0.0,0.0,0.0
	write(9,"('H',/,3F20.14)") ax(1),ax(2),ax(3)
	write(9,"('O',/,3F20.14)") ay(1),ay(2),ay(3)
	write(9,"('S',/,3F20.14)") az(1),az(2),az(3)
	end program ortho

	subroutine getvector(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	implicit none
	real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3
	! xyz1,xyz2: coordinates of two atoms; xyz3: vector between them (result)
	x3=x1-x2
	y3=y1-y2
	z3=z1-z2
	end subroutine getvector
