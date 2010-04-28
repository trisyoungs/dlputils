	! **dlx2** program to convert between formatted and unformatted DL_POLY hisfiles
	EXTERNAL readframe,writeframe,readheader,writeheader
	integer maxn, readframe, writeframe,readheader,writeheader
	PARAMETER(MAXN=8000)
	character*80 header,infile,outfile
	character*8 atmnam(MAXN)
	character*8 discard
	integer success, outform
	integer keytrj, imcon, natms, nstep, iatm, n, discardn, nframes, l, fileform
	real*8 tstep
	real*8 cell(9), xpos(MAXN),ypos(MAXN),zpos(MAXN), charge(MAXN), mass(MAXN)
	real*8 xvel(MAXN),yvel(MAXN),zvel(MAXN),xfor(MAXN),yfor(MAXN),zfor(MAXN)
	real*8 volume, avvolume, maxvolume,minvolume, mina,maxa,minb,maxb,minc,maxc
	real*8 ava,avb,avc
	COMMON /frame/ atmnam,xpos,ypos,zpos,charge,mass,xvel,yvel,zvel,xfor,yfor,zfor,
     +  cell,header,tstep,keytrj,imcon,natms,nstep


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
	  outfile=infile(1:l-1)//"f"
	  OPEN(UNIT=11,FILE=outfile,FORM="FORMATTED")
	ELSE
	  OPEN(UNIT=10,FILE=infile,FORM="FORMATTED")
	  outfile=infile(1:l-1)//"u"
	  OPEN(UNIT=11,FILE=outfile,FORM="UNFORMATTED")
	END IF

	! Do the conversion.
	! Read/write header...
	IF (readheader(10,fileform).EQ.-1) GOTO 979
	IF (writeheader(11,outform).EQ.-1) GOTO 978
100     success=readframe(10,fileform)
	IF (success.EQ.1) GOTO 980  ! End of file encountered....
	IF (success.EQ.-1) GOTO 979  ! File error....
	IF (writeframe(11,outform).EQ.-1) GOTO 978
	IF (MOD(nstep,10000).LT.10) WRITE(0,*) infile(1:l),": ",nstep
	GOTO 100
	

978     WRITE(0,*) "Error while writing to output file (step",nstep,")."
	GOTO 999
979	WRITE(0,*) "History file ended prematurely at step ",nstep,"!"
	GOTO 999
980	WRITE(*,*) "Finished."
999	CLOSE(10)
	CLOSE(11)
	END
	
