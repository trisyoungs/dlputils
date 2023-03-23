	! Change the atom names in a CONFIG file to the users spec.
	! Can do batch processing of a number of files.

	program editconfig
	use dlprw
	implicit none
	integer :: edittype,editatom,o,p,success,nfiles,n,m,count,f
	integer :: iargc
        character*80 :: outfile
	character*80, allocatable :: filenames(:)
	character*8, allocatable :: oldnames(:,:),newnames(:,:)
	character*8 :: t1,t2

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

	! Get the list of target CONFIG files...
	nfiles = iargc()
	if (nfiles.EQ.0) stop "No files specified!"
	allocate(filenames(nfiles))
	do n=1,nfiles
	  call getarg(n,filenames(n))
	end do

	! Open the OUTPUT file if it exists, or get user input.
	write(6,*) "Name of OUTPUT file:"
	read(5,*) outfile
	success = outinfo(outfile,1)

	allocate(oldnames(nspecies,maxatoms))
	allocate(newnames(nspecies,maxatoms))

	! Read in the 'old' atom names from the first CONFIG file....
	call readconfig(filenames(1))
	count = 0
	do n=1,nspecies
	 do o=1,s_nmols(n)
	  do m=1,s_natoms(n)
	    count = count + 1
	    oldnames(n,m)=cfgname(count)
	    newnames(n,m)=oldnames(n,m)
	  end do
	 end do
	end do

	! Print the list of atom names now and ask for changes......
70	write(*,*) ""
	write(*,*) "Select molecular species to edit (99 to exit):"
	read(*,*) edittype
	if (edittype.EQ.99) goto 80
	if (edittype.gt.nspecies) goto 70
71	write(*,*) "Atom names of molecular species ",edittype,":"
	write(*,*) "    Original    New"
	do n=1,s_natoms(edittype)
	  if (oldnames(edittype,n).EQ.newnames(edittype,n)) THEN
	    write(*,"(I2,4X,A5)") n,oldnames(edittype,n)
	  else
	    write(*,"(I2,4X,A5,A5,A5)") n,oldnames(edittype,n)," --> ",newnames(edittype,n)
	  end if
	end do
	write(*,*) ""
	write(*,*) "Enter atom number to change (or zero for a different species):"
	read(*,*) editatom
	if (editatom.EQ.0) goto 70 
	if (editatom.LT.0) THEN
	  write(*,*) "Change all atoms of type:"
	  read(*,"(A5)") t1
	  write(*,*) "To?"
	  read(*,"(A5)") t2
	  do n=1,s_natoms(edittype)
	    if (oldnames(edittype,n).EQ.t1) newnames(edittype,n)=t2
	  end do
	else
	  write(*,*) "New name:"
	  read(*,*) newnames(edittype,editatom)
	end if
	write(*,*) ""
	goto 71

	! Main loop begins here.
	! Have the list of changes to make, so read frame, make changes, write frame...
80	write(6,"(A,I3,A)") "Processing ",nfiles," CONFIGs..."
	do f=1,nfiles
	  count = 0
	  write(6,"(A,A)") "Modifying ",filenames(f)
	  call readconfig(filenames(f))
	  do n=1,nspecies
	    do m=1,s_nmols(n)
	      do o=1,s_natoms(n)
	        count = count + 1
	        cfgname(count) = newnames(n,o)
	      end do
	    end do
	  end do
	  call writeconfig(filenames(f))
	  write(6,"(A,A)") "Finished."
	end do

	goto 990

	! Exit messages.
200	write(*,*) "Error while reading in CONFIG file data."
	write(*,*) "Were the molecular species specified correctly?"
	do p=1,nspecies
	  write(*,*) "  Type ",p,": nmols=",s_nmols(p),"  natoms=",s_natoms(p)
	end do
	write(*,*) "Error occured on Type(",n,") Mol(",o,") Atom(",m,")."
990	stop
	end program editconfig
