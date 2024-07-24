	! Routines to probe DL_POLY unformatted history files.

	program probedlp
	implicit none
	! History file
	character*80 :: hisfile
	character*8 discard, atmname
	integer :: n, success, nargs
	! History file data
	character*80 :: dlpheader, temp
	real*8 :: cell(9), tstep, x, x1, x2, x3, x4, x5
	integer :: keytrj, imcon, natms, nstep
	integer :: iargc

	! Get target filename
	nargs = iargc()
	if (nargs.ne.1) stop "Usage: probedlp <HISTORYfile>"
	call getarg(1,hisfile)

	write(0,"(A,A)") "History file : ",hisfile

	! Open history file
	open(unit=40,file=hisfile,form='unformatted',status='old')

	! Attempt to read header
	success = 0
	write(0,*) "Reading header from unformatted history file..."
	! File dlpheader
	read(40,ERR=91,end=91) dlpheader
	goto 92
91	write(0,*) "Failed to read 80-character header string."
	goto 200

	! Number of atoms, followed by individual names, masses and charges
92	read(40,ERR=93,end=93) x
	natms=IDNINT(x)
	write(0,*) "Number of atoms:",natms
	goto 94
93	write(0,*) "Failed to read number of atoms from header."
	goto 200
94	read(40,ERR=95,end=95) (atmname,n=1,natms)
	write(0,*) "Read atomnames successfully."
	goto 96
95	write(0,*) "Failed to read atom name", n
	goto 200

96	read(40,ERR=97,end=97) (x,n=1,natms)
	write(0,*) "Read atom masses successfully."
	goto 98
97	write(0,*) "Failed to read atom mass", n
	goto 200

98	read(40,ERR=99,end=99) (x,n=1,natms)
	write(0,*) "Read atom charges successfully."
	success = 1
	goto 200
99	write(0,*) "Failed to read atom charge", n
	goto 200

200	if (success.eq.0) then
	  write(0,*) "Failed to read header info. Rewinding to start of file."
	  rewind(40)
	end if

	! Now ready for frame data

300	success = 0

	! Step number, natoms, keytrj, imcon, and timestep 
	x1 = 0.0
	x2 = 0.0
	x3 = 0.0
	x4 = 0.0
	x5 = 0.0
	read(40,ERR=301,end=500) x1,x2,x3,x4,x5
	n=IDNINT(x1)
	keytrj=IDNINT(x3)
	imcon=IDNINT(x4)
	write(0,*) "Step number:",n
	write(0,*) "Number of atoms:",IDNINT(x2)
	write(0,*) "Trajectory key:",keytrj
	write(0,*) "Image Convention:",imcon
	write(0,*) "Timestep:",x5
	goto 302
301	write(0,*) "Failed to read step/natms/keytrj/imcon/timestep from frame."
	write(0,*) x1,x2,x3,x4,x5
	goto 400

	! Cell
302	if (imcon.gt.0) then
	  read(40,ERR=303,end=303) (x,n=1,9)
	  write(0,*) "Read cell matrix."
	  goto 304
303	  write(0,*) "Failed to read cell parameter", n
	  goto 400
	end if
	   
	! X-Coordinates
304	read(40,ERR=305,end=305) (x,n=1,natms)
	write(0,*) "Read x-coordinates."
	goto 306
305	write(0,*) "Failed to read x-coordinate", n
	goto 400

	! Y-Coordinates
306	read(40,ERR=307,end=307) (x,n=1,natms)
	write(0,*) "Read y-coordinates."
	goto 308
307	write(0,*) "Failed to read y-coordinate", n
	goto 400

	! Z-Coordinates
308	read(40,ERR=309,end=309) (x,n=1,natms)
	write(0,*) "Read z-coordinates."
	goto 310
309	write(0,*) "Failed to read z-coordinate", n
	goto 400

310	if (keytrj.gt.0) then
	  ! X-Velocities
	  read(40,ERR=311,end=311) (x,n=1,natms)
	  write(0,*) "Read x-velocities."
	  goto 312
311	  write(0,*) "Failed to read x-velocity", n
	  goto 400

	  ! Y-velocities
312	  read(40,ERR=313,end=313) (x,n=1,natms)
	  write(0,*) "Read y-velocities."
	  goto 314
313	  write(0,*) "Failed to read y-velocity", n
	  goto 400

	  ! Z-velocities
314	  read(40,ERR=315,end=315) (x,n=1,natms)
	  write(0,*) "Read z-velocities."
	  goto 316
315	  write(0,*) "Failed to read z-velocity", n
	  goto 400
	end if

316	if (keytrj.gt.1) then
	  ! X-forces
	  read(40,ERR=317,end=317) (x,n=1,natms)
	  write(0,*) "Read x-forces."
	  goto 318
317	  write(0,*) "Failed to read x-velocity", n
	  goto 400

	  ! Y-forces
318	  read(40,ERR=319,end=319) (x,n=1,natms)
	  write(0,*) "Read y-forces."
	  goto 320
319	  write(0,*) "Failed to read y-velocity", n
	  goto 400

	  ! Z-forces
320	  read(40,ERR=321,end=321) (x,n=1,natms)
	  write(0,*) "Read z-forces."
	  goto 380 
321	  write(0,*) "Failed to read z-velocity", n
	end if

380	success = 1
400	if (success.eq.1) goto 300
	write(0,*) "Failed to read frame."
	stop
500	write(0,*) "Found logical end of file."

	end program probedlp
