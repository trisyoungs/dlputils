	! **getsysdata** extracts system data from *.OUT and writes a file
	implicit none
	! Variables
	character*300 discard,dis
	character*29 discard2
	integer step,n,m,numsteps,interval,numpass,oddpass,numdata
	real eng_tot,temp_tot,eng_cfg,eng_vdw,eng_cou,eng_bnd,eng_ang
	real eng_dih,eng_tet,eng_pv,press,discardr,volume,temp_shl,eng_shl
	
10	FORMAT (9X,9(1X,E11.4E2))
11	FORMAT (I9,9(1X,E11.4E2))
14	FORMAT (I9,11(1X,E11.4E2))
15	FORMAT (12(A10,2X))

	WRITE(6,15) "Step","E-Total","Temp  ","E-config","E-vdw  ","E-cou  ","E-bond ","E-angle","E-torsion","E-PV  ","Pressure", "E-shell"
	! Now we start the main loop
100	READ(5,11,END=200,ERR=200) step,eng_tot,temp_tot,eng_cfg,eng_vdw,eng_cou,eng_bnd,eng_ang,eng_dih,eng_tet
	READ(5,10,END=200,ERR=200) eng_pv,discardr,discardr,discardr,discardr,discardr,discardr,discardr,discardr
	READ(5,10,END=200,ERR=200) volume,temp_shl,eng_shl,discardr,discardr,discardr,discardr,discardr,press
	WRITE(6,14) step,eng_tot,temp_tot,eng_cfg,eng_vdw,eng_cou,eng_bnd,eng_ang,eng_dih,eng_pv,press,eng_shl
	GOTO 100

200	WRITE(0,*) "Complete."
	CLOSE(5)
	CLOSE(6)
	END
