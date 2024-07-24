!	** glucgeom **
!	Analyses the interactions in our simulations

	program glucgeom
	use dlprw; use utility
	implicit none
	character*80 :: hisfile,dlpoutfile,resfile, flagfile
	character*20 :: temparg
	character*4 :: cutstr
	integer :: n,nframes,success,nargs,m,i,mol,glucoff,cloff
	integer :: framestodo, nsites, t(9), nglucose, nchloride, gluc, cl, nframesused, gluc2
	integer, allocatable :: sites(:,:), cl_ohsite(:,:), molflags(:), glucosecn(:), siteflags(:,:)
	integer, allocatable :: glucgluc11(:,:), glucgluc21(:,:), glucgluc22(:,:), ggmap(:), hb_per_gluc(:)
	integer :: hist_distbins, hist_anglebins, startf=1, m1, m2, getpair, m3, m4
	integer :: iargc, ggmappairs(3)
	logical :: molmap = .FALSE., writeggmap = .FALSE.
	real*8 :: total, count_cn, count_cloh(3), count_hb, count_gg11, count_gg21, count_gg22
	real*8 :: a(3),b(3),c(3),temp(3),v1(3),v2(3),v1mag,v2mag,roo,angle
	real*8, allocatable :: cl_ohdist(:,:), cl_ohangle(:,:)
	real*8 :: calcangle, cutoff
	real*8 :: hist_distbin, hist_anglebin, hist_distmax, hist_anglemax
	real*8, allocatable :: dist_cloh(:,:), angle_cloh(:,:), angle_hclh_intra(:), angle_hclh_inter(:), site_cl(:,:)
	real*8, allocatable :: site_ohpairs_inter(:,:), site_ohpairs_intra(:,:)

	nargs = iargc()
	if (nargs.LT.3) stop "Usage : glucgeom <HISTORYfile> <OUTfile> <cutoff> [-mapfile file] [-start frame] [-frames n] [-writeggmap n m 1|2]"
	call getarg(1,hisfile)
	call getarg(2,dlpoutfile)
	call getarg(3,temparg); read(temparg,"(f10.5)") cutoff
	n = 3
	do
	  n = n + 1; if (n.GT.nargs) exit
	  call getarg(n,temparg)
	  select case (temparg)
	    case ("-frames")
	      n = n + 1; call getarg(n,temparg); read(temparg,"(I6)") framestodo
	      write(0,"(A,I6)") "Number of frames to sample: ",framestodo
	    case ("-start")
	      n = n + 1; call getarg(n,temparg); read(temparg,"(I6)") startf
	      write(0,"(A,I4)") "Starting frame = ",startf
	    case ("-molmap")
	      n = n + 1; call getarg(n,flagfile)
	      write(0,"(A,A)") "Using flags for species 1 molecules from : ",flagfile
	      open(unit=25,file=flagfile,form="formatted",status="old")
	      molmap = .TRUE.
	    case ("-writeggmap")
	      n = n + 1; call getarg(n,temparg); read(temparg,"(i4)") m
	      n = n + 1; call getarg(n,temparg); read(temparg,"(i4)") i
	      n = n + 1; call getarg(n,temparg); read(temparg,"(i4)") ggmappairs(3)
	      ggmappairs(1) = min(m,i)
	      ggmappairs(2) = max(m,i)
	      write(0,"(A,i2,'-',i2)") "Will write mapfile for gluc=Cl=gluc interactions."
	      write(0,"(A,i2,'-',i2,a,i2)") " -- Pairs are ",ggmappairs(1:2),", centering on ",ggmappairs(3)
	      writeggmap = .TRUE.
	    case default
	      write(0,*) "Unrecognised command line option : ",temparg
	      stop
	  end select
	end do
	
	! Make string of cutoff distance
	write(cutstr,"(f4.1)") cutoff
	if (int(cutoff).lt.10) cutstr(1:1) = "0"

	! Open and check the files...
	call openhis(hisfile,10)
	if (outinfo(dlpoutfile,1).EQ.-1) goto 798
	! Now, read in the history header so that we have cell()
	if (readheader().EQ.-1) goto 799
	! Open the glucose CN and sitemap files ready for writing
	open(unit=30,file="glucgeom-glucCN-"//cutstr//".dat",form="formatted",status="replace")
	open(unit=31,file="glucgeom-sitemap-"//cutstr//".dat",form="formatted",status="replace")
	if (writeggmap) open(unit=32,file="glucgeom-glucglucmap.dat",form="formatted",status="replace")

	! Numbers
	nsites = 5
	nglucose = s_nmols(1)
	nchloride = s_nmols(3)
	hist_distbin = 0.05
	hist_distmax = 15.0
	hist_anglebin = 0.5
	hist_anglemax = 180.0
	hist_distbins = hist_distmax / hist_distbin
	hist_anglebins = hist_anglemax / hist_anglebin

	! Data
	allocate(molflags(nglucose))
	allocate(sites(nsites,2))
	allocate(cl_ohdist(nglucose,nsites))
	allocate(cl_ohangle(nglucose,nsites))
	allocate(cl_ohsite(nglucose,nsites))

	! Analysis
	allocate(siteflags(nglucose,nsites))
	allocate(hb_per_gluc(nglucose))
	allocate(glucosecn(nglucose))
	allocate(dist_cloh(3,hist_distbins))
	allocate(angle_cloh(3,hist_anglebins))
	allocate(angle_hclh_intra(hist_anglebins))
	allocate(angle_hclh_inter(hist_anglebins))
	allocate(site_cl(2,nsites))
	allocate(site_ohpairs_inter(nsites,nsites))
	allocate(site_ohpairs_intra(nsites,nsites))
	allocate(glucgluc11(nsites,nsites))
	allocate(glucgluc21(nsites*nsites,nsites))
	allocate(glucgluc22(nsites*nsites,nsites*nsites))
	allocate(ggmap(nglucose))

	! Initialise
	dist_cloh = 0.0
	angle_cloh = 0.0
	angle_hclh_intra = 0.0
	angle_hclh_inter = 0.0
	site_cl = 0
	site_ohpairs_inter = 0
	site_ohpairs_intra = 0
	count_cloh = 0
	count_cn = 0
	count_hb = 0
	count_gg11 = 0
	count_gg21 = 0
	count_gg22 = 0
	total = 0
	glucgluc11 = 0
	glucgluc21 = 0
	glucgluc22 = 0

	!sites(:,:) = reshape( (/ 8,9,10,11,12,20,21,22,23,24 /) , (/ nsites,2 /) )
	sites(:,:) = reshape( (/ 9,10,11,12,8,21,22,23,24,20 /) , (/ nsites,2 /) )

	write(0,*) "OH site specification (O, H):"
	do m=1,nsites
	  write(0,*) (sites(m,n),n=1,2)
	end do
	write(0,"(A,A)") "Name of chloride species (check) : ",s_name(3)

	! XXXX
	! XXXX Main routine....
	! XXXX
	! Set up the vars...
100	nframes=0
	nframesused = 0
101     if (molmap) then
	  read(25,"(I5,1x,40(I2,1x))",end=118,err=117) m1,(molflags(m2),m2=1,nglucose)
	  ! if (m1.NE.nframes) stop "Mapfile missed a frame number!"
	end if
102	success=readframe()
	if (success.EQ.1) goto 801  ! End of file encountered....
	if (success.LT.0) goto 799  ! File error....
	nframes=nframes+1
        if (mod(nframes,100).EQ.0) then
          if (molmap) then
            write(0,"(I6,2x,'(',I6,',',f8.0,')')") nframes,nframesused,total
          else
            write(0,"(i6)") nframes
          end if
        end if
	if (nframes.lt.startf) goto 102
	if ((molmap).and.(m1.ne.nframes)) goto 102

	! ANALYSIS BEGINS HERE
	glucosecn = 0
	siteflags = 0
	ggmap = 0
	hb_per_gluc = 0
	nframesused = nframesused + 1

	! Loop over chlorides in the outer loop so we can calc all Cl-HO geometries in one go
	cloff = s_start(3) - 1
	do cl=1,nchloride

	  ! Grab chloride position
	  c(1)=xpos(cloff+cl)
	  c(2)=ypos(cloff+cl)
	  c(3)=zpos(cloff+cl)
	  ! Calculate all Cl-HO interactions for this chloride
	  glucoff = 0
	  do gluc=1,nglucose
	    do m=1,nsites
		a(1)=xpos(sites(m,1)+glucoff)	! Glucose Oxygen
		a(2)=ypos(sites(m,1)+glucoff)
		a(3)=zpos(sites(m,1)+glucoff)
		b(1)=xpos(sites(m,2)+glucoff)	! Glucose Hydrogen
		b(2)=ypos(sites(m,2)+glucoff)
		b(3)=zpos(sites(m,2)+glucoff)
		! Cl-O distance
		call pbc(a(1),a(2),a(3),c(1),c(2),c(3),temp(1),temp(2),temp(3))
		cl_ohdist(gluc,m) = sqrt( (temp(1)-c(1))**2 + (temp(2)-c(2))**2 + (temp(3)-c(3))**2 )
		if (cl_ohdist(gluc,m).lt.2.0) write(0,*) "Warning!!! Short Chloride-OH contact.",sites(m,1),cloff+n,nframes
		! Cl-H-O angle
		cl_ohangle(gluc,m) = calcangle(a,b,c)
		cl_ohsite(gluc,m) = m
	    end do	! Loop over glucose OH sites
	    ! Sort by distance
	    call sort(nsites,cl_ohdist(gluc,:),cl_ohangle(gluc,:),cl_ohsite(gluc,:))
	    glucoff = glucoff + 24
	  end do	! Loop over glucose molecules

	  ! Now we have all geometries we need to do the analysis per glucose molecule
	  do gluc=1,nglucose
	    if (molmap.and.(molflags(gluc).eq.0)) cycle
	    total = total + 1	! General counter

	    ! Is the current chloride within 'cutoff' angstroms of one of the OH groups of the current glucose?
	    ! Check shortest Cl-HO distance to find out
	    if (cl_ohdist(gluc,1).le.cutoff) then

	      ! Increase global and glucose-specific CN counters
	      count_cn = count_cn + 1
	      glucosecn(gluc) = glucosecn(gluc) + 1
	      ! Bin the Cl-HO angle for this, the closest Cl-H interaction of this particular chloride
	      m = (cl_ohangle(gluc,1) / hist_anglebin) + 1
	      if (m.le.hist_anglebins) angle_cloh(3,m) = angle_cloh(3,m) + 1
	
	      ! Is it also within 'cutoff' distance from a second OH on the same molecule?
	      if (cl_ohdist(gluc,2).le.cutoff) then
		! Three?!?
		if (cl_ohdist(gluc,3).le.cutoff) then
		  ! Could start to panic here....
		  count_cloh(3) = count_cloh(3) + 1
		else
		  ! An intramolecular bridging chloride (presumably)
		  count_cloh(2) = count_cloh(2) + 1
		  ! Add distances and angles to histogram
		  m = (cl_ohdist(gluc,1) / hist_distbin) + 1
		  if (m.le.hist_distbins) dist_cloh(2,m) = dist_cloh(2,m) + 1
		  m = (cl_ohdist(gluc,2) / hist_distbin) + 1
		  if (m.le.hist_distbins) dist_cloh(2,m) = dist_cloh(2,m) + 1
		  m = (cl_ohangle(gluc,1) / hist_anglebin) + 1
		  if (m.le.hist_anglebins) angle_cloh(2,m) = angle_cloh(2,m) + 1
		  m = (cl_ohangle(gluc,2) / hist_anglebin) + 1
		  if (m.le.hist_anglebins) angle_cloh(2,m) = angle_cloh(2,m) + 1
		  ! Store OH groups involved and get H-Cl-H coordinates to calc angle
		  glucoff = (gluc-1) * 24
		  m1 = cl_ohsite(gluc,1)
		  site_cl(2,m1) = site_cl(2,m1) + 1
		  siteflags(gluc,m1) = siteflags(gluc,m1) + 1
	  	  a(1)=xpos(sites(m1,2)+glucoff)
	  	  a(2)=ypos(sites(m1,2)+glucoff)
	  	  a(3)=zpos(sites(m1,2)+glucoff)
		  m2 = cl_ohsite(gluc,2)
		  site_cl(2,m2) = site_cl(2,m2) + 1
		  siteflags(gluc,m2) = siteflags(gluc,m2) + 1
	  	  b(1)=xpos(sites(m2,2)+glucoff)
	  	  b(2)=ypos(sites(m2,2)+glucoff)
	  	  b(3)=zpos(sites(m2,2)+glucoff)
		  site_ohpairs_intra(m1,m2) = site_ohpairs_intra(m1,m2) + 1
		  site_ohpairs_intra(m2,m1) = site_ohpairs_intra(m2,m1) + 1
		  ! Get H-Cl-H angle of interaction and add to a histogram
	  	  c(1)=xpos(s_start(3)+cl-1)
	  	  c(2)=ypos(s_start(3)+cl-1)
	  	  c(3)=zpos(s_start(3)+cl-1)
		  angle = calcangle(a,c,b)
		  m = (angle / hist_anglebin) + 1
		  angle_hclh_intra(m) = angle_hclh_intra(m) + 1
		  ! Here, check whether another glucose is bound to the chloride
		  do gluc2=gluc+1,nglucose
		    if (cl_ohdist(gluc2,1).le.cutoff) then
		      ! There is a chloride involved in a bridge with 'gluc', and at least one H-bond with gluc2...
		      ! So, check if there is a bridge on gluc2...
		      if (cl_ohdist(gluc2,2).le.cutoff) then
			! The chloride is involved in four HBonds, two with each gluc...
			! Find the id of the *pair* of OHs that are involved for each glucose
			m1 = getpair(cl_ohsite(gluc,1),cl_ohsite(gluc,2))
			m2 = getpair(cl_ohsite(gluc2,1),cl_ohsite(gluc2,2))
			m3 = min(m1,m2)
			m4 = max(m1,m2)
			glucgluc22(m3,m4) = glucgluc22(m3,m4) + 1
			count_gg22 = count_gg22 + 1
			if ((m3.eq.ggmappairs(1)).and.(m4.eq.ggmappairs(2))) then
			  ! Pair matches, so find which one we want to flag
			  if (ggmappairs(3).eq.1) then
			    ggmap(gluc) = 1
			  else
			    ggmap(gluc2) = 1
			  end if
			end if
		      else
			! The chloride is involved in a bridge with 'gluc' OHs and s single HBond with 'gluc2'
			m1 = getpair(cl_ohsite(gluc,1),cl_ohsite(gluc,2))
			m2 = cl_ohsite(gluc2,1)
			glucgluc21(m1,m2) = glucgluc21(m1,m2) + 1
			count_gg21 = count_gg21 + 1
			write(89,*) nframes
		      end if
		    end if
		  end do
		end if
	      else
		! No, only one close OH contact (to this glucose)
		count_cloh(1) = count_cloh(1) + 1
		! Add distance and angle to histogram
		m = (cl_ohdist(gluc,1) / hist_distbin) + 1
		if (m.le.hist_distbins) dist_cloh(1,m) = dist_cloh(1,m) + 1
		m = (cl_ohangle(gluc,1) / hist_anglebin) + 1
		if (m.le.hist_anglebins) angle_cloh(1,m) = angle_cloh(1,m) + 1
		! Store OH group involved
		m = cl_ohsite(gluc,1)
		site_cl(1,m) = site_cl(1,m) + 1
		siteflags(gluc,m) = siteflags(gluc,m) + 2

		! Check if short contact exists to other glucose molecules
		do gluc2=gluc+1,nglucose
		  ! Only need to check for single OH interaction - 2-1 interactions taken care of in upper loop
		  if (cl_ohdist(gluc2,1).le.cutoff) then
		    m3 = cl_ohsite(gluc,1)
		    m4 = cl_ohsite(gluc2,1)
		    m1 = min(m3,m4)
		    m2 = max(m3,m4)
		    glucgluc11(m1,m2) = glucgluc11(m1,m2) + 1
		    count_gg11 = count_gg11 + 1
		    write(90,*) nframes
		  end if
		end do

	      end if	! Second interaction cutoff test
	    end if	! First interaction cutoff test

	    ! Independent HB total
	    do m=1,nsites
	      if (cl_ohdist(gluc,m).le.cutoff) then
		count_hb = count_hb + 1
		hb_per_gluc(gluc) = hb_per_gluc(gluc) + 1
	      end if
	    end do

	  end do	! Loop over glucose molecules

	end do	! Loop over chlorides

	! Write out the glucosecn and sitemaps for use in other things
	write(30,"(i5,2x,40I2)") nframes,glucosecn
	write(31,"(i5,2x,20(5i1,1x))") nframes,siteflags
	if (writeggmap) write(32,"(i5,1x,40(i2,1x))") nframes,ggmap
	write(44,"(i2)") hb_per_gluc

	if (nframesused.EQ.framestodo) goto 800
	! Next frame
	goto 101

117     write(0,*) "Error reading mapfile."
        write(0,*) "Selected ",nframesused," from trajectory before error."
        goto 801
118     write(0,*) "Reached end of mapfile."
        write(0,*) "Selected ",nframesused," from trajectory."
        goto 801
700	write(*,*) "INFILE and OUTFILE have the same name!!!"
	goto 999
798	write(0,*) "Problem with OUTPUT file."
	goto 999
799	write(0,*) "HISTORY file ended before framestodo was fulfilled..."
	write(0,"(A,I4,A,I4,A)") "Frames read in : ",nframesused," (wanted ",framestodo,")"
	goto 801
800	write(0,*) "Framestodo was fulfilled."
801	write(0,*) ""

	total = total / nchloride

	write(0,"(a,f7.2)") "	              Distance cutoff : ",cutoff
	if (molmap) then
	  write(0,"(a,a)") "	                      Mapfile : ",flagfile
	else
	  write(0,"(a)") "	                      Mapfile : None"
	end if
	write(0,"(a,i5)") "	                  Start frame : ",startf
	write(0,"(a,i5)") "	            Total frames used : ",nframesused
	write(0,"(a,f7.0)") "Total glucose environments considered : ",total
	write(0,*) ""
	write(0,"(a,f7.4)") "   Average CN (per glucose per frame) : ",real(count_cn / total)
	write(0,"(a,f7.0,' (',f5.1,'%)')") 		"	 Number of Cl bound to one OH : ",count_cloh(1),(count_cloh(1) / count_cn)*100.0
	write(0,"(a,5(f7.0,1x),' (',5(f5.1,1x),'%)')")  "	             Site populations : ",(site_cl(1,m),m=1,5),((site_cl(1,m)/count_cloh(1))*100.0,m=1,5)
	write(0,"(a,f7.0,' (',f5.1,'%)')") "	 Number of Cl bound to two OH : ",count_cloh(2),(count_cloh(2) / count_cn)*100.0
	write(0,"(a,5(f7.0,1x),' (',5(f5.1,1x),'%)')")  "	             Site populations : ",(site_cl(2,m),m=1,5),((site_cl(2,m)/count_cloh(2))*50.0,m=1,5)
	write(0,"(a,f7.0,' (',f5.1,'%)')") "       Number of Cl bound to three OH : ",count_cloh(3),(count_cloh(3) / count_cn)*100.0
	!write(0,"(a,5(f7.0,1x),' (',5(f5.1,1x),'%)')")  "	             Site populations : ",(site_cl(3,m),m=1,5),((site_cl(3,m)/count_cloh(3))*100.0,m=1,5)
	write(0,"(a,f7.4)") " HBond total (from above populations) : ",(count_cloh(1)+count_cloh(2)*2+count_cloh(3)*3) / total
	write(0,"(a,f7.4)") "    Independent HBond total (<cutoff) : ",count_hb / total

	open(unit=14,file="glucgeom-"//cutstr//".out",form="formatted",status="replace")
	write(14,"(a,f7.2)") "	              Distance cutoff : ",cutoff
	if (molmap) then
	  write(14,"(a,a)") "	                      Mapfile : ",flagfile
	else
	  write(14,"(a)") "	                      Mapfile : None"
	end if
	write(14,"(a,i5)") "	                  Start frame : ",startf
	write(14,"(a,i5)") "	            Total frames used : ",nframesused
	write(14,"(a,f7.0)") "Total glucose environments considered : ",total
	write(14,*) ""
	write(14,"(a,f7.4)") "   Average CN (per glucose per frame) : ",real(count_cn / total)
	write(14,"(a,f7.0,' (',f5.1,'%)')") "	 Number of Cl bound to one OH : ",count_cloh(1),(count_cloh(1) / count_cn)*100.0
	write(14,"(a,5(f7.0,1x),' (',5(f5.1,1x),'%)')")  "	             Site populations : ",(site_cl(1,m),m=1,5),((site_cl(1,m)/count_cloh(1))*100.0,m=1,5)
	write(14,"(a,f7.0,' (',f5.1,'%)')") "	 Number of Cl bound to two OH : ",count_cloh(2),(count_cloh(2) / count_cn)*100.0
	write(14,"(a,5(f7.0,1x),' (',5(f5.1,1x),'%)')")  "	             Site populations : ",(site_cl(2,m),m=1,5),((site_cl(2,m)/count_cloh(2))*50.0,m=1,5)
	write(14,"(a,f7.0,' (',f5.1,'%)')") "	Number of Cl bound to three OH : ",count_cloh(3),(count_cloh(3) / count_cn)*100.0
	!write(14,"(a,5(f7.0,1x),' (',5(f5.1,1x),'%)')")  "	             Site populations : ",(site_cl(3,m),m=1,5),((site_cl(3,m)/count_cloh(3))*100.0,m=1,5)
	write(14,"(a,f7.4)") " HBond total (from above populations) : ",(count_cloh(1)+count_cloh(2)*2+count_cloh(3)*3) / total
	write(14,"(a,f7.4)") "    Independent HBond total (<cutoff) : ",count_hb / total
	write(14,*) ""
	write(14,"(a)") "    Glucose-chloride-glucose interactions:"
	write(14,"(a,f7.0,' (',f7.4,')')") "         Number of OH-Cl-HO interactions (p2gpf) : ",count_gg11,2*count_gg11 / total
	write(14,"(a,f7.0,' (',f7.4,')')") "        Number of Bridge-HO interactions (p2gpf) : ",count_gg21,2*count_gg21 / total
	write(14,"(a,f7.0,' (',f7.4,')')") "    Number of Bridge-Bridge interactions (p2gpf) : ",count_gg22,2*count_gg22 / total

	! Write out histograms etc.

	open(unit=14,file="glucgeom-disthist-"//cutstr//".dat",form="formatted",status="replace")
	write(14,"(a)") "#   r      cl-1oh    cl-2oh    cl-3oh  "
	do n=1,hist_distbins
	  write(14,"(f7.4,2x,3f10.2)") n*hist_distbin-hist_distbin/2.0, (dist_cloh(m,n),m=1,3)
	end do
	close(14)
	
	open(unit=14,file="glucgeom-anglehist-"//cutstr//".dat",form="formatted",status="replace")
	write(14,"(a)") "# angle    cl-1oh    cl-2oh    closest   hclh(intra)  hclh(inter)"
	do n=1,hist_anglebins
	  write(14,"(f7.3,2x,4f10.2)") n*hist_anglebin-hist_anglebin/2.0,(angle_cloh(m,n),m=1,3),angle_hclh_intra(n)
	end do
	close(14)
	
	open(unit=14,file="glucgeom-sitepops-"//cutstr//".dat",form="formatted",status="replace")
	write(14,"(a)") "# site    oh-cl     oh-cl-oh    oh-cl     oh-cl-oh     "
	do n=1,5
	  write(14,"(i5,2x,2f10.0,2f10.4)") n,(site_cl(m,n),m=1,2),(site_cl(m,n)/count_cloh(m)*(100/m),m=1,2)
	end do
	close(14)
	
	open(unit=14,file="glucgeom-ohpairs-intra-"//cutstr//".dat",form="formatted",status="replace")
	write(14,"(4x,5i10)") (m,m=1,5)
	do n=1,5
	  write(14,"(i3,2x,5f10.4)") n,(site_ohpairs_intra(n,m)/count_cloh(2)*100.0,m=1,5)
	end do
	close(14)

	open(unit=14,file="glucgluc-bridgebridge-"//cutstr//".dat",form="formatted",status="replace")
	write(14,"(5x,25(1x,i1,'-',i1,2x))") ((m-1)/5+1,mod(m-1,5)+1,m=1,25)
	do n=1,25
	  write(14,"(i1,'-',i1,2x,25(f5.2,1x))") (n-1)/5+1,mod(n-1,5)+1,(glucgluc22(n,m)/count_gg22*100.0,m=1,25)
	end do
	close(14)

	open(unit=14,file="glucgluc-bridgeoh-"//cutstr//".dat",form="formatted",status="replace")
	write(14,"(4x,5i10)") (m,m=1,5)
	do n=1,25
	  write(14,"(i1,'-',i1,2x,5f10.4)") (n-1)/5+1,mod(n-1,5)+1,(glucgluc21(n,m)/count_gg21*100.0,m=1,5)
	end do
	close(14)

	open(unit=14,file="glucgluc-ohoh-"//cutstr//".dat",form="formatted",status="replace")
	write(14,"(4x,5i10)") (m,m=1,5)
	do n=1,5
	  write(14,"(i3,2x,5f10.4)") n,(glucgluc11(n,m)/count_gg11*100.0,m=1,5)
	end do
	close(14)

	write(0,*) "Finished."
999	close(14)
	end program glucgeom

	integer function getpair(i,j)
	implicit none
	integer :: i,j,ii,jj
	! Create unique ID for pair of OHs (i and j) that are involved 
	! site order = 8,9,10,11,12
	ii = min(i,j)
	jj = max(i,j)
	getpair = (ii-1)*5 + jj
	!getpair = 5	!Default
	!if ((ii.eq.1).and.(jj.eq.5)) getpair = 4
	!if ((ii.eq.2).and.(jj.eq.3)) getpair = 1
	!if ((ii.eq.3).and.(jj.eq.4)) getpair = 2
	!if ((ii.eq.4).and.(jj.eq.5)) getpair = 3
	end function getpair
	
	real*8 function calcangle(a,b,c)
	use utility; implicit none
	real*8,dimension(3) :: a,b,c,v1,v2,temp1,temp2
	real*8 :: v1mag,v2mag,dp
	  ! Angle between atoms a--b--c
	call pbc(a(1),a(2),a(3),b(1),b(2),b(3),temp1(1),temp1(2),temp1(3))
	call pbc(c(1),c(2),c(3),b(1),b(2),b(3),temp2(1),temp2(2),temp2(3))
	! Get the bond vectors
	v1 = b - temp1
	v2 = b - temp2
	v1mag = sqrt( v1(1)**2 + v1(2)**2 + v1(3)**2 )
	v2mag = sqrt( v2(1)**2 + v2(2)**2 + v2(3)**2 )
	v1 = v1 / v1mag
	v2 = v2 / v2mag
	  ! Calculate dot product and angle...
	  dp = (v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3))
	  calcangle = safeAngle(dp)
	end function calcangle


	SUBROUTINE sort(n,arr,brr,crr)
	implicit none
	INTEGER n,M,NSTACK
	DOUBLE PRECISION arr(n),brr(n)
	PARAMETER (M=7,NSTACK=50)
	INTEGER i,ir,j,jstack,k,l,istack(NSTACK),c,crr(n)
	DOUBLE PRECISION a,b,temp,tempi
	jstack=0
	l=1
	ir=n
1     if(ir-l.lt.M)then
	  do 12 j=l+1,ir
	    a=arr(j)
	    b=brr(j)
	    c=crr(j)
	    do 11 i=j-1,1,-1
	      if(arr(i).le.a)goto 2
	      arr(i+1)=arr(i)
	      brr(i+1)=brr(i)
	      crr(i+1)=crr(i)
11	  continue
	    i=0
2	   arr(i+1)=a
	    brr(i+1)=b
	    crr(i+1)=c
12	continue
	  if(jstack.eq.0)return
	  ir=istack(jstack)
	  l=istack(jstack-1)
	  jstack=jstack-2
	else
	  k=(l+ir)/2
	  temp=arr(k)
	  arr(k)=arr(l+1)
	  arr(l+1)=temp
	  temp=brr(k)
	  brr(k)=brr(l+1)
	  brr(l+1)=temp
	  tempi=crr(k)
	  crr(k)=crr(l+1)
	  crr(l+1)=tempi
	  if(arr(l+1).gt.arr(ir))then
	    temp=arr(l+1)
	    arr(l+1)=arr(ir)
	    arr(ir)=temp
	    temp=brr(l+1)
	    brr(l+1)=brr(ir)
	    brr(ir)=temp
	    tempi=crr(l+1)
	    crr(l+1)=crr(ir)
	    crr(ir)=tempi
	  endif
	  if(arr(l).gt.arr(ir))then
	    temp=arr(l)
	    arr(l)=arr(ir)
	    arr(ir)=temp
	    temp=brr(l)
	    brr(l)=brr(ir)
	    brr(ir)=temp
	    tempi=crr(l)
	    crr(l)=crr(ir)
	    crr(ir)=tempi
	  endif
	  if(arr(l+1).gt.arr(l))then
	    temp=arr(l+1)
	    arr(l+1)=arr(l)
	    arr(l)=temp
	    temp=brr(l+1)
	    brr(l+1)=brr(l)
	    brr(l)=temp
	    tempi=crr(l+1)
	    crr(l+1)=crr(l)
	    crr(l)=tempi
	  endif
	  i=l+1
	  j=ir
	  a=arr(l)
	  b=brr(l)
	  c=crr(l)
3	 continue
	    i=i+1
	  if(arr(i).lt.a)goto 3
4	 continue
	    j=j-1
	  if(arr(j).gt.a)goto 4
	  if(j.lt.i)goto 5
	  temp=arr(i)
	  arr(i)=arr(j)
	  arr(j)=temp
	  temp=brr(i)
	  brr(i)=brr(j)
	  brr(j)=temp
	  tempi=crr(i)
	  crr(i)=crr(j)
	  crr(j)=tempi
	  goto 3
5	 arr(l)=arr(j)
	  arr(j)=a
	  brr(l)=brr(j)
	  brr(j)=b
	  crr(l)=crr(j)
	  crr(j)=c
	  jstack=jstack+2
	  if(jstack.gt.NSTACK) stop "Sort2 error"
	  if(ir-i+1.ge.j-l)then
	    istack(jstack)=ir
	    istack(jstack-1)=i
	    ir=j-1
	  else
	    istack(jstack)=j-1
	    istack(jstack-1)=l
	    l=i
	  endif
	endif
	goto 1
	END
