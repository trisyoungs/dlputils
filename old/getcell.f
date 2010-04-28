!	** getcell **
!	Retrieves the cell dimensions from a HISTORY file and outputs them to a file suitable
!	for plotting in Excel, along with a summary and averages.
	EXTERNAL readframe,readheader
	integer maxn,readframe,readheader
	PARAMETER(MAXN=8000)
	character*80 header,infile,outfile
	character*8 atmnam(MAXN)
	character*8 discard
	integer success
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
	ELSE IF (infile(l:l).EQ."F") THEN
	  fileform=1
	ELSE IF (infile(l:l).EQ."u") THEN
	  fileform=0
	ELSE IF (infile(l:l).EQ."U") THEN
	  fileform=0
	ELSE
	  WRITE(*,*) "Unable to determine the form of the HISTORY file."
	  WRITE(*,*) "Please enter 0 for unformatted or 1 for formatted:"
	  READ(*,*) fileform
	END IF

	IF (fileform.EQ.0) THEN
	  OPEN(UNIT=10,FILE=infile,FORM="UNFORMATTED")
	ELSE
	  OPEN(UNIT=10,FILE=infile,FORM="FORMATTED")
	END IF

	WRITE(*,*) "Name of output file?"
	READ(*,*) outfile
	IF (infile.EQ.outfile) GOTO 700
	OPEN(UNIT=11,FILE=outfile,FORM="FORMATTED")

	! Have what we need, so read in the stuff from before to get to the start of the trajectory.....
	IF (readheader(10,fileform).EQ.-1) GOTO 799
	! Also, set the necessary variables to zero.....
	avvolume=0
	minvolume=1E20
	maxvolume=0
	mina=1E5
	minb=1E5
	minc=1E5
	maxa=0
	maxb=0
	maxc=0
	ava=0
	avb=0
	avc=0
	nframes=0
	WRITE(11,"(A58)") "    Step          a           b           c       Volume   "

	! At START of trajectory
100	success=readframe(10,fileform)
	IF (success.EQ.1) GOTO 800  ! End of file encountered....
	IF (success.EQ.-1) GOTO 799  ! File error....
	! Calculate the current box volume and update averages e.t.c.
	volume=cell(1)*cell(5)*cell(9)
	avvolume=avvolume+volume/1E4
	if (volume.LT.minvolume) minvolume=volume
	if (volume.GT.maxvolume) maxvolume=volume
	ava=ava+cell(1)
	avb=avb+cell(5)
	avc=avc+cell(9)
	if (cell(1).LT.mina) mina=cell(1)
	if (cell(1).GT.maxa) maxa=cell(1)
	if (cell(5).LT.minb) minb=cell(5)
	if (cell(5).GT.maxb) maxb=cell(5)
	if (cell(9).LT.minc) minc=cell(9)
	if (cell(9).GT.maxc) maxc=cell(9)
	! Now write out the data.....
	WRITE(11,"(I8,2X,4F12.4)") nstep,cell(1),cell(5),cell(9),volume
	nframes=nframes+1
	GOTO 100

700	WRITE(*,*) "INFILE and OUTFILE have the same name!!!"
	GOTO 999

799	WRITE(11,"HISTORY file ended prematurely!")
	GOTO 801
800	WRITE(11,*) "End of unformatted HISTORY file found."
801	WRITE(11,*) "Final averages......"
	ava=ava/nframes
	avb=avb/nframes
	avc=avc/nframes
	avvolume=(avvolume/nframes)*1E4
	WRITE(11,"(A58)") "                  a           b           c       Volume   "
	WRITE(11,"(A10,4F12.4)") "Minimum",mina,minb,minc,minvolume
	WRITE(11,"(A10,4F12.4)") "Maximum",maxa,maxb,maxc,maxvolume
	WRITE(11,"(A10,4F12.4)") "Average",ava,avb,avc,avvolume
	
999	CLOSE(unitno)
	CLOSE(11)
	END

