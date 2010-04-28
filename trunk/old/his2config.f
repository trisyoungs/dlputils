	! ** his2config ** reads a history frame and writes a config file.
	EXTERNAL readframe,writeframe,readheader,writeheader
	integer maxn, readframe, writeframe,readheader,writeheader
	PARAMETER(MAXN=8000)
	character*80 header,infile,outfile
	character*8 atmnam(MAXN)
	character*8 discard
	integer success, outform, target
	integer keytrj, imcon, natms, nstep, iatm, n, discardn, nframes, l, fileform
	real*8 tstep
	real*8 cell(9), xpos(MAXN),ypos(MAXN),zpos(MAXN), charge(MAXN), mass(MAXN)
	real*8 xvel(MAXN),yvel(MAXN),zvel(MAXN),xfor(MAXN),yfor(MAXN),zfor(MAXN)
	real*8 volume, avvolume, maxvolume,minvolume, mina,maxa,minb,maxb,minc,maxc
	real*8 ava,avb,avc
	COMMON /frame/ atmnam,xpos,ypos,zpos,charge,mass,xvel,yvel,zvel,xfor,yfor,zfor,
     +  cell,header,tstep,keytrj,imcon,natms,nstep

15      FORMAT (A5/,3f20.14)
16      FORMAT (3(2x,F18.10))

	WRITE(*,*) "Name of HISTORY file?"
	READ(*,*) infile
	! Determine the FORM of the HISTORY file...
	l=INDEX(infile," ")-1
	IF (infile(l:l).EQ."f") THEN
	  fileform=1
	  outform=0
	ELSE IF (infile(l:l).EQ."F") THEN
	  fileform=1
	  outform=0
	ELSE IF (infile(l:l).EQ."u") THEN
	  fileform=0
	  outform=1
	ELSE IF (infile(l:l).EQ."U") THEN
	  fileform=0
	  outform=1
	ELSE
	  WRITE(*,*) "Unable to determine the form of the HISTORY file."
	  WRITE(*,*) "Please enter 0 for unformatted or 1 for formatted:"
	  READ(*,*) fileform
	  IF(fileform.EQ.1) THEN
	    outform=0
	  ELSE
	    outform=1
	  END IF
	END IF

	IF (fileform.EQ.0) THEN
	  OPEN(UNIT=10,FILE=infile,FORM="UNFORMATTED")
	ELSE
	  OPEN(UNIT=10,FILE=infile,FORM="FORMATTED")
	END IF

	! Get the target frame....
	WRITE(0,*) "Which frame to convert?"
	READ(*,*) target

	IF (readheader(10,fileform).EQ.-1) GOTO 979
	DO n=1,target
	  success=readframe(10,fileform)
	  IF (success.EQ.1) GOTO 979  ! End of file encountered....
	  IF (success.EQ.-1) GOTO 979  ! File error....
	END DO
	
	! Write out the config file....
	WRITE(0,*) imcon, keytrj, natms
	outfile="frame.CONFIG"
	OPEN(UNIT=13,FILE=outfile,FORM="FORMATTED")
	WRITE(13,*) "Frame output from his2config."
	IF (imcon.GT.0) THEN
	  WRITE(13,"(2I10)") keytrj,imcon
	  WRITE(13,16) cell(1),0.0,0.0,0.0,cell(5),0.0,0.0,0.0,cell(9)
	ELSE
	  WRITE(13,"(2I10)") keytrj,0
	END IF
	DO n=1,natms
	  WRITE(13,15) atmnam(n),xpos(n),ypos(n),zpos(n)
	  WRITE(0,15) atmnam(n),xpos(n),ypos(n),zpos(n)
	  IF (keytrj.GT.0) WRITE(13,16) xvel(n),yvel(n),zvel(n)
	  IF (keytrj.GT.1) WRITE(13,16) xfor(n),yfor(n),zfor(n)
	END DO
	CLOSE(13)
	GOTO 999


978     WRITE(0,*) "Error while writing to output file (step",nstep,")."
	GOTO 999
979	WRITE(0,*) "History file ended prematurely at step ",nstep,"!"
	GOTO 999
980	WRITE(*,*) "Finished."
999	CLOSE(10)
	CLOSE(13)
	END
	
