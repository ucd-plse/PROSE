!******************************************************************************
!  MODULE RS2
!    Written by s.b. 09/xx/2006
!------------------------------------------------------------------------------
! This module reads in UNIT 23 file and external files specified in UNIT 23.
! This module was made to have ADCIRC read in a single global UNIT 23 file
! and sets of global STWAVE output files (.sim and .rad files) from each 
! processor.
!------------------------------------------------------------------------------
! sb  12/11/2006  RS2GET now takes maximum values if stwave domains overlap.
!                 (It used to take the value of the domain of the first 
!                 appearance in a fort.23 file.)
! vjp 10/21/2006  modified logic for type2
!****************************************************************************** 
! 
      MODULE RS2
      USE SIZES

      USE GLOBAL,ONLY : NSCREEN, ScreenUnit
#ifdef CMPI
      use messenger, only : msg_fini
#endif
      IMPLICIT NONE

      !variables for input file
      integer :: ns,type
      character(len=100),allocatable,dimension(:) :: simfiles ! names of .sim files
      character(len=100),allocatable,dimension(:) :: radfiles ! names of .rad files
      integer,allocatable,dimension(:)            :: fullplanes ! 0: half plane, 1: full plane
      integer,allocatable,dimension(:)            :: nbs ! number of blank snaps
      integer,allocatable,dimension(:)            :: rfids ! id of .rad files
      integer,allocatable,dimension(:) :: nis,njs
      integer,save,allocatable,dimension(:) :: endoffile

      integer,allocatable,dimension(:,:,:) :: rsij
      real(kind=8),allocatable,dimension(:,:,:) :: rscnf

      real(kind=8),allocatable,dimension(:) :: dxincs,dyincs
      real(kind=8),allocatable,dimension(:) :: x0s,y0s,azimuths

      integer :: numSets,numBlankSnaps,cntSnaps,numSkipSnaps
      REAL(kind=8) :: windMultiplier


      PUBLIC

!---------------------end of data declarations--------------------------------C


      CONTAINS


!***********************************************************************
!   SOBROUTINE RS2INIT
!***********************************************************************

      subroutine rs2init(rsnx,rsny,np)
      use sizes
      implicit none
      integer,intent(in) :: np
      real(kind=8),intent(out) :: rsnx(:),rsny(:)

      if(myproc.eq.0) then
        write(screenunit,*)
        write(screenunit,*) 'INIT RADIATION STRESS ARRAYS (NRS=2)...'
      endif
      write(16,*)
      write(16,*) 'INITIALIZING RADIATION STRESS ARRAYS (NRS=2)...'

      !read meta info


      OPEN(23,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.23',STATUS='OLD')

      ! type of fort.23 file, 1: include wave radiation
      !                       2: does not include wave radiation
      read(23,*,err=99999) type

      if(type.eq.1) then
        write(*,'(a)') 'RS2INIT: SORRY! TYPE-1 RS FILE NOT YET SUPPORTED.'
        CALL EXIT(1)
!        call rs2init_type1(rsnx,rsny,np)
      else if(type.eq.2) then
        call rs2init_type2(rsnx,rsny,np)
      else
         write(*,'(a)') 'RS2 TYPE IS NEITHER 1 NOR 2.  TERMINATED.'
         CALL EXIT(1)
      endif

      return

99999 CONTINUE

#ifdef CMPI
      call msg_fini()
#endif
      WRITE(*,'(a)') 'RADIATION STRESS READ ERROR (RS2-0)'
      CALL EXIT(1)

      END SUBROUTINE


!$$$C***********************************************************************
!$$$C   SOBROUTINE RS2INIT_TYPE1
!$$$C***********************************************************************
!$$$
!$$$      subroutine rs2init_type1(rsnx,rsny,np)
!$$$      use sizes
!$$$      use global,only : slam,sfea,rad2deg
!$$$      implicit none
!$$$      integer,intent(in) :: np
!$$$      real(kind=8),intent(out) :: rsnx(:),rsny(:)
!$$$      integer ::  c,i,j,k,n,p,s,ni,nj
!$$$      real(kind=8) :: lon,lat
!$$$      real(kind=8) :: dxinc,dyinc
!$$$      real(kind=8) :: stwlon1,stwlat1,stwlon2,stwlat2
!$$$      real(kind=8) :: stwlon3,stwlat3,stwlon4,stwlat4
!$$$      real(kind=8) :: x1,y1,x2,y2,x3,y3
!$$$      real(kind=8) :: subarea1,subarea2,subarea3,subarea4,totalarea
!$$$      real(kind=8),dimension(:,:,:), allocatable :: stwlonlat
!$$$      real(kind=8) :: stwlonmin,stwlonmax,stwlatmin,stwlatmax
!$$$      integer,parameter :: nSearchBins = 50
!$$$      integer,allocatable :: nCellsInSearchBins(:)
!$$$      integer,allocatable :: searchBins(:,:,:)
!$$$      real(kind=8),allocatable :: searchBinPartitions(:)
!$$$      real(kind=8) :: p1,p2,cmax,cmin,tol
!$$$
!$$$
!$$$      write(*,'(a)') 'RS2INIT_TYPE1: SORRY! TYPE-1 RS FILE NOT YET SUPPORTED.'
!$$$      CALL EXIT(1) 
!$$$
!$$$      read(23,*,err=99999) ns    !number of stwave grid files
!$$$
!$$$      if(ns.le.0) then
!$$$        write(screenunit,1004)
!$$$        write(16,1004)
!$$$#ifdef CMPI
!$$$        call msg_fini()
!$$$#endif
!$$$        CALL EXIT(1)
!$$$      endif
!$$$
!$$$      allocate(fullplanes(ns),rfids(ns))
!$$$      allocate(nis(ns),njs(ns),dxincs(ns),dyincs(ns))
!$$$      allocate(x0s(ns),y0s(ns),azimuths(ns))
!$$$      allocate(endoffile(ns))
!$$$
!$$$      allocate(rsij(4,np),rscnf(3,np))
!$$$
!$$$      do n=1,np
!$$$        rsij(1,n) = 0
!$$$      enddo
!$$$
!$$$      !compute mapping coefficients
!$$$      do s=1,ns
!$$$        read(23,*,err=99999) ! skip a line
!$$$        read(23,*,err=99999) fullplanes(s)
!$$$
!$$$        if(fullplanes(s).eq.1) then
!$$$          read(23,*,err=99999) nis(s), njs(s), dxincs(s), dyincs(s)
!$$$        else
!$$$          read(23,*,err=99999) nis(s), njs(s), dxincs(s)
!$$$          dyincs(s) = dxincs(s)
!$$$        endif       
!$$$        
!$$$        ni = nis(s)
!$$$        nj = njs(s)
!$$$
!$$$        allocate(stwlonlat(2,ni,nj))
!$$$        do j=1,nj
!$$$          read(23,*) ((stwlonlat(k,i,j), k=1,2), i=1,ni)
!$$$        enddo
!$$$        
!$$$        stwlonmin =  1e10
!$$$        stwlonmax = -1e10
!$$$        stwlatmin =  1e10
!$$$        stwlatmax = -1e10
!$$$        do j=1,nj
!$$$          do i=1,ni
!$$$            if(stwlonlat(1,i,j).eq.0.e0) cycle
!$$$            if(stwlonmin.gt.stwlonlat(1,i,j)) stwlonmin=stwlonlat(1,i,j)
!$$$            if(stwlonmax.lt.stwlonlat(1,i,j)) stwlonmax=stwlonlat(1,i,j)
!$$$            if(stwlatmin.gt.stwlonlat(2,i,j)) stwlatmin=stwlonlat(2,i,j)
!$$$            if(stwlatmax.lt.stwlonlat(2,i,j)) stwlatmax=stwlonlat(2,i,j)
!$$$          enddo
!$$$        enddo
!$$$
!$$$        !prepare search bins
!$$$        allocate(nCellsInSearchBins(nSearchBins))
!$$$        allocate(searchBins(2,(ni-1)*(nj-1),nSearchBins))
!$$$        allocate(searchBinPartitions(nSearchBins+1))
!$$$        searchBinPartitions(1) = stwlonmin
!$$$        if(myproc.eq.0) then
!$$$          write(screenunit,*) 'PREPARING SEARCH BINS FOR GRID ',s,'...'
!$$$        endif
!$$$        tol = abs(stwlonmax-stwlonmin)/real(nSearchBins)*0.01e0
!$$$        do p=1,nSearchBins
!$$$          searchBinPartitions(p+1) = 
!$$$     &         stwlonmin + (stwlonmax-stwlonmin)*
!$$$     &         real(p)/real(nSearchBins)
!$$$
!$$$          p1 = searchBinPartitions(p)
!$$$          p2 = searchBinPartitions(p+1)
!$$$
!$$$          p1 = p1 - tol
!$$$          p2 = p2 + tol
!$$$
!$$$          nCellsInSearchBins(p) = 0
!$$$
!$$$          do j=1,nj-1
!$$$            do i=1,ni-1
!$$$              stwlon1 = stwlonlat(1,i,j)
!$$$              stwlon2 = stwlonlat(1,i+1,j)
!$$$              stwlon3 = stwlonlat(1,i+1,j+1)
!$$$              stwlon4 = stwlonlat(1,i,j+1)
!$$$
!$$$              cmax = stwlon1
!$$$              if(cmax.lt.stwlon2) cmax = stwlon2
!$$$              if(cmax.lt.stwlon3) cmax = stwlon3
!$$$              if(cmax.lt.stwlon4) cmax = stwlon4
!$$$
!$$$              cmin = stwlon1
!$$$              if(cmin.gt.stwlon2) cmin = stwlon2
!$$$              if(cmin.gt.stwlon3) cmin = stwlon3
!$$$              if(cmin.gt.stwlon4) cmin = stwlon4
!$$$
!$$$              if(cmax.ge.p1.and.cmin.le.p2) then
!$$$                 nCellsInSearchBins(p) = nCellsInSearchBins(p) + 1
!$$$                 searchBins(1,nCellsInSearchBins(p),p) = i
!$$$                 searchBins(2,nCellsInSearchBins(p),p) = j
!$$$              endif
!$$$            enddo
!$$$          enddo
!$$$        enddo
!$$$        if(myproc.eq.0) then
!$$$          write(screenunit,*) ' SEARCH BINS ARE READY'
!$$$        endif
!$$$
!$$$        do n=1,np
!$$$
!$$$          if(rsij(1,n).ne.0) cycle
!$$$
!$$$          lat = RAD2DEG*SFEA(n)
!$$$          lon = RAD2DEG*SLAM(n)
!$$$
!$$$          if(lon.lt.stwlonmin.or.lon.gt.stwlonmax.or.
!$$$     &       lat.lt.stwlatmin.or.lat.gt.stwlatmax) cycle
!$$$            
!$$$
!$$$          do p=1,nSearchBins
!$$$            p1 = searchBinPartitions(p)
!$$$            p2 = searchBinPartitions(p+1)
!$$$            
!$$$            if(lon.ge.p1.and.lon.le.p2) exit
!$$$          enddo
!$$$
!$$$          if(p.gt.nSearchBins) cycle
!$$$
!$$$          cellloop: do c=1,nCellsInSearchBins(p)
!$$$            i = searchBins(1,c,p)
!$$$            j = searchBins(2,c,p)
!$$$            stwlon1 = stwlonlat(1,i,j)
!$$$            stwlat1 = stwlonlat(2,i,j)
!$$$            stwlon2 = stwlonlat(1,i+1,j)
!$$$            stwlat2 = stwlonlat(2,i+1,j)
!$$$            stwlon3 = stwlonlat(1,i+1,j+1)
!$$$            stwlat3 = stwlonlat(2,i+1,j+1)
!$$$            stwlon4 = stwlonlat(1,i,j+1)
!$$$            stwlat4 = stwlonlat(2,i,j+1)
!$$$
!$$$            if(stwlon1.eq.0.e0.or.stwlon2.eq.0.e0.or.
!$$$     &           stwlon3.eq.0.e0.or.stwlon4.eq.0.e0) then
!$$$              cycle
!$$$            endif
!$$$
!$$$            !triangle 1 (nodes 1, 2 and 3)
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon2
!$$$            y2 = stwlat2
!$$$            x3 = stwlon3
!$$$            y3 = stwlat3
!$$$            subarea1 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$            
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon3
!$$$            y2 = stwlat3
!$$$            x3 = stwlon1
!$$$            y3 = stwlat1
!$$$            subarea2 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$            
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon1
!$$$            y2 = stwlat1
!$$$            x3 = stwlon2
!$$$            y3 = stwlat2
!$$$            subarea3 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$            
!$$$            x1 = stwlon1
!$$$            y1 = stwlat1
!$$$            x2 = stwlon2
!$$$            y2 = stwlat2
!$$$            x3 = stwlon3
!$$$            y3 = stwlat3
!$$$            totalarea = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$            
!$$$            if((subarea1+subarea2+subarea3).le.
!$$$     &           (totalarea*1.00001e0)) then
!$$$              rsij(1,n) = s
!$$$              rsij(2,n) = 1
!$$$              rsij(3,n) = i
!$$$              rsij(4,n) = j
!$$$              rscnf(1,n)=subarea1/totalarea
!$$$              rscnf(2,n)=subarea2/totalarea
!$$$              rscnf(3,n)=subarea3/totalarea
!$$$              exit cellloop
!$$$            endif
!$$$
!$$$            !triangle 2 (nodes 3, 4 and 1)
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon4
!$$$            y2 = stwlat4
!$$$            x3 = stwlon1
!$$$            y3 = stwlat1
!$$$            subarea1 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$            
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon1
!$$$            y2 = stwlat1
!$$$            x3 = stwlon3
!$$$            y3 = stwlat3
!$$$            subarea2 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$            
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon3
!$$$            y2 = stwlat3
!$$$            x3 = stwlon4
!$$$            y3 = stwlat4
!$$$            subarea3 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$            
!$$$            x1 = stwlon3
!$$$            y1 = stwlat3
!$$$            x2 = stwlon4
!$$$            y2 = stwlat4
!$$$            x3 = stwlon1
!$$$            y3 = stwlat1
!$$$            totalarea = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$            
!$$$            if((subarea1+subarea2+subarea3).le.
!$$$     &           (totalarea*1.00001e0)) then
!$$$              rsij(1,n) = s
!$$$              rsij(2,n) = 2
!$$$              rsij(3,n) = i
!$$$              rsij(4,n) = j
!$$$              rscnf(1,n)=subarea1/totalarea
!$$$              rscnf(2,n)=subarea2/totalarea
!$$$              rscnf(3,n)=subarea3/totalarea
!$$$              exit cellloop
!$$$            endif
!$$$          enddo cellloop
!$$$        enddo
!$$$
!$$$        deallocate(stwlonlat)
!$$$        deallocate(nCellsInSearchBins)
!$$$        deallocate(searchBins)
!$$$        deallocate(searchBinPartitions)
!$$$
!$$$      enddo
!$$$
!$$$      !initialize some variables
!$$$      cntSnaps = 0
!$$$      endoffile(:) = 0
!$$$
!$$$      if(myproc.eq.0) then
!$$$        write(screenunit,*)
!$$$        write(screenunit,*) 'INIT OF RADIATION STRESS ARRAYS IS DONE.'
!$$$      endif
!$$$      write(16,*)
!$$$      write(16,*) 'INIT OF RADIATION STRESS ARRAYS IS DONE.'
!$$$
!$$$      return
!$$$
!$$$ 1004 FORMAT(//,1X,' NUMBER OF SETS WAS SPECIFIED'//
!$$$     &       'INCORRECTLY IN UNIT23.'/
!$$$     &       ' IT MUST BE GREATER THAN 0.'//)
!$$$
!$$$99999 CONTINUE
!$$$
!$$$#ifdef CMPI
!$$$      call msg_fini()
!$$$#endif
!$$$      write(*,'(a)') 'RADIATION STRESS READ ERROR (RS2-1)'
!$$$      CALL EXIT(1) 
!$$$
!$$$      END SUBROUTINE


!***********************************************************************
!   SOBROUTINE RS2INIT_TYPE2
!***********************************************************************

      subroutine rs2init_type2(rsnx,rsny,np)
      use sizes
      use global,only : rad2deg
      use mesh, only : slam, sfea
      implicit none
      integer,intent(in) :: np
      real(kind=8),intent(out) :: rsnx(:),rsny(:)
      integer ::  c,i,j,k,n,p,s,ni,nj
      real(kind=8) :: lon,lat
      real(kind=8) :: dxinc,dyinc
      real(kind=8) :: stwlon1,stwlat1,stwlon2,stwlat2
      real(kind=8) :: stwlon3,stwlat3,stwlon4,stwlat4
      real(kind=8) :: x1,y1,x2,y2,x3,y3
      real(kind=8) :: subarea1,subarea2,subarea3,subarea4,totalarea
      real(kind=8),dimension(:,:,:), allocatable :: stwlonlat
      real(kind=8) :: stwlonmin,stwlonmax,stwlatmin,stwlatmax
      integer,parameter :: nSearchBins = 50
      integer,allocatable :: nCellsInSearchBins(:)
      integer,allocatable :: searchBins(:,:,:)
      real(kind=8),allocatable :: searchBinPartitions(:)
      real(kind=8) :: p1,p2,cmax,cmin,tol


      read(23,*,err=99999) ns    !number of stwave grid files

      if(ns.le.0) then
        write(screenunit,1004)
        write(16,1004)
#ifdef CMPI
        call msg_fini()
#endif
        CALL EXIT(1)
      endif

      allocate(fullplanes(ns),rfids(ns))
      allocate(nis(ns),njs(ns),dxincs(ns),dyincs(ns))
      allocate(x0s(ns),y0s(ns),azimuths(ns))
      allocate(endoffile(ns))
      allocate(simfiles(ns),radfiles(ns))
      allocate(nbs(ns))

      allocate(rsij(3,np,ns),rscnf(3,np,ns))

      rsij(1,:,:) = 0
         

      !compute mapping coefficients
      do s=1,ns
        read(23,*,err=99999) ! skip a line

        read(23,'(A)',err=99999) simfiles(s)
        simfiles(s) = adjustl(simfiles(s))
        read(23,'(A)',err=99999) radfiles(s)
        radfiles(s) = adjustl(radfiles(s))
        read(23,*,err=99999) nbs(s)
        if (myproc == 0) then
         write(screenunit,*) " NUMBER BLANK RS SNAPS FOR REGION ",s,
     &                         " = ", nbs(s)
        endif
        read(23,*,err=99999) x0s(s),y0s(s),azimuths(s)

        read(23,*,err=99999) fullplanes(s)

        if(fullplanes(s).eq.1) then
          read(23,*,err=99999) nis(s), njs(s), dxincs(s), dyincs(s)
        else
          read(23,*,err=99999) nis(s), njs(s), dxincs(s)
          dyincs(s) = dxincs(s)
        endif       
        
        ni = nis(s)
        nj = njs(s)

        allocate(stwlonlat(2,ni,nj))
        do j=1,nj
          read(23,*) ((stwlonlat(k,i,j), k=1,2), i=1,ni)
        enddo
        
        stwlonmin =  1e10
        stwlonmax = -1e10
        stwlatmin =  1e10
        stwlatmax = -1e10
        do j=1,nj
          do i=1,ni
            if(stwlonlat(1,i,j).eq.0.e0) cycle
            if(stwlonmin.gt.stwlonlat(1,i,j)) stwlonmin=stwlonlat(1,i,j)
            if(stwlonmax.lt.stwlonlat(1,i,j)) stwlonmax=stwlonlat(1,i,j)
            if(stwlatmin.gt.stwlonlat(2,i,j)) stwlatmin=stwlonlat(2,i,j)
            if(stwlatmax.lt.stwlonlat(2,i,j)) stwlatmax=stwlonlat(2,i,j)
          enddo
        enddo

        !prepare search bins
        allocate(nCellsInSearchBins(nSearchBins))
        allocate(searchBins(2,(ni-1)*(nj-1),nSearchBins))
        allocate(searchBinPartitions(nSearchBins+1))
        searchBinPartitions(1) = stwlonmin
        if(myproc.eq.0) then
          write(screenunit,*) 'PREPARING SEARCH BINS FOR GRID ',s,'...'
        endif
        tol = abs(stwlonmax-stwlonmin)/real(nSearchBins)*0.01e0
        do p=1,nSearchBins
          searchBinPartitions(p+1) = 
     &         stwlonmin + (stwlonmax-stwlonmin)*
     &         real(p)/real(nSearchBins)

          p1 = searchBinPartitions(p)
          p2 = searchBinPartitions(p+1)

          p1 = p1 - tol
          p2 = p2 + tol

          nCellsInSearchBins(p) = 0

          do j=1,nj-1
            do i=1,ni-1
              stwlon1 = stwlonlat(1,i,j)
              stwlon2 = stwlonlat(1,i+1,j)
              stwlon3 = stwlonlat(1,i+1,j+1)
              stwlon4 = stwlonlat(1,i,j+1)

              cmax = stwlon1
              if(cmax.lt.stwlon2) cmax = stwlon2
              if(cmax.lt.stwlon3) cmax = stwlon3
              if(cmax.lt.stwlon4) cmax = stwlon4

              cmin = stwlon1
              if(cmin.gt.stwlon2) cmin = stwlon2
              if(cmin.gt.stwlon3) cmin = stwlon3
              if(cmin.gt.stwlon4) cmin = stwlon4

              if(cmax.ge.p1.and.cmin.le.p2) then
                 nCellsInSearchBins(p) = nCellsInSearchBins(p) + 1
                 searchBins(1,nCellsInSearchBins(p),p) = i
                 searchBins(2,nCellsInSearchBins(p),p) = j
              endif
            enddo
          enddo
        enddo
        if(myproc.eq.0) then
          write(screenunit,*) ' SEARCH BINS ARE READY'
        endif

        do n=1,np

          lat = RAD2DEG*SFEA(n)
          lon = RAD2DEG*SLAM(n)

          if(lon.lt.stwlonmin.or.lon.gt.stwlonmax.or.
     &       lat.lt.stwlatmin.or.lat.gt.stwlatmax) cycle
            

          do p=1,nSearchBins
            p1 = searchBinPartitions(p)
            p2 = searchBinPartitions(p+1)
            
            if(lon.ge.p1.and.lon.le.p2) exit
          enddo

          if(p.gt.nSearchBins) cycle

          cellloop: do c=1,nCellsInSearchBins(p)
            i = searchBins(1,c,p)
            j = searchBins(2,c,p)
            stwlon1 = stwlonlat(1,i,j)
            stwlat1 = stwlonlat(2,i,j)
            stwlon2 = stwlonlat(1,i+1,j)
            stwlat2 = stwlonlat(2,i+1,j)
            stwlon3 = stwlonlat(1,i+1,j+1)
            stwlat3 = stwlonlat(2,i+1,j+1)
            stwlon4 = stwlonlat(1,i,j+1)
            stwlat4 = stwlonlat(2,i,j+1)

            if(stwlon1.eq.0.e0.or.stwlon2.eq.0.e0.or.
     &           stwlon3.eq.0.e0.or.stwlon4.eq.0.e0) then
              cycle
            endif

            !triangle 1 (nodes 1, 2 and 3)
            x1 = lon
            y1 = lat
            x2 = stwlon2
            y2 = stwlat2
            x3 = stwlon3
            y3 = stwlat3
            subarea1 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
            
            x1 = lon
            y1 = lat
            x2 = stwlon3
            y2 = stwlat3
            x3 = stwlon1
            y3 = stwlat1
            subarea2 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
            
            x1 = lon
            y1 = lat
            x2 = stwlon1
            y2 = stwlat1
            x3 = stwlon2
            y3 = stwlat2
            subarea3 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
            
            x1 = stwlon1
            y1 = stwlat1
            x2 = stwlon2
            y2 = stwlat2
            x3 = stwlon3
            y3 = stwlat3
            totalarea = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
            
            if((subarea1+subarea2+subarea3).le.
     &           (totalarea*1.00001e0)) then
              rsij(1,n,s) = 1
              rsij(2,n,s) = i
              rsij(3,n,s) = j
              rscnf(1,n,s)=subarea1/totalarea
              rscnf(2,n,s)=subarea2/totalarea
              rscnf(3,n,s)=subarea3/totalarea
              exit cellloop
            endif

            !triangle 2 (nodes 3, 4 and 1)
            x1 = lon
            y1 = lat
            x2 = stwlon4
            y2 = stwlat4
            x3 = stwlon1
            y3 = stwlat1
            subarea1 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
            
            x1 = lon
            y1 = lat
            x2 = stwlon1
            y2 = stwlat1
            x3 = stwlon3
            y3 = stwlat3
            subarea2 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
            
            x1 = lon
            y1 = lat
            x2 = stwlon3
            y2 = stwlat3
            x3 = stwlon4
            y3 = stwlat4
            subarea3 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
            
            x1 = stwlon3
            y1 = stwlat3
            x2 = stwlon4
            y2 = stwlat4
            x3 = stwlon1
            y3 = stwlat1
            totalarea = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
            
            if((subarea1+subarea2+subarea3).le.
     &           (totalarea*1.00001e0)) then
              rsij(1,n,s) = 2
              rsij(2,n,s) = i
              rsij(3,n,s) = j
              rscnf(1,n,s)=subarea1/totalarea
              rscnf(2,n,s)=subarea2/totalarea
              rscnf(3,n,s)=subarea3/totalarea
              exit cellloop
            endif
          enddo cellloop
        enddo

        deallocate(stwlonlat)
        deallocate(nCellsInSearchBins)
        deallocate(searchBins)
        deallocate(searchBinPartitions)

      enddo

      !open .rad files if necessary
      do s=1,ns
        rfids(s) = 230 + s
        open(rfids(s),file=TRIM(GBLINPUTDIR)//'/'//trim(radfiles(s)),
     & status='old')
        if(fullplanes(s).eq.1) then
          read (rfids(s),*) ni, nj, dxinc, dyinc
        else
          read (rfids(s),*) ni, nj, dxinc
          dyinc = dxinc
        endif
        if((ni.ne.nis(s)).or.(nj.ne.njs(s)).or.
     &       (dxinc.ne.dxincs(s)).or.(dyinc.ne.dyincs(s))) then
          if(myproc.eq.0) then
            write(screenunit,*) ''
            write(screenunit,*) 'FORT.23 INCONSISTENT WITH .RAD FILES'
            write(screenunit,*) 'PROGRAM WILL BE TERMINATED'
            write(screenunit,*) ''
          endif
          write(16,*) ''
          write(16,*) 'FORT.23 SEEMS INCONSISTENT WITH .RAD FILES'
          write(16,*) 'PROGRAM WILL BE TERMINATED'
          write(16,*) ''
#ifdef CMPI
          call msg_fini()
#endif
          CALL EXIT(1)
        endif
      enddo

      !skip snaps if nbs(s) < 0
      do s=1,ns
        do
          if(nbs(s).ge.0) exit
            nbs(s) = nbs(s) + 1
            call rs2get(RSNX,RSNY,NP)
        enddo
      enddo

      !initialize some variables
      cntSnaps = 0
      endoffile(:) = 0

      if(myproc.eq.0) then
        write(screenunit,*)
        write(screenunit,*) 'INIT OF RADIATION STRESS ARRAYS IS DONE.'
      endif
      write(16,*)
      write(16,*) 'INIT OF RADIATION STRESS ARRAYS IS DONE.'

      return

 1004 FORMAT(//,1X,' NUMBER OF SETS WAS SPECIFIED'//
     &       'INCORRECTLY IN UNIT23.'/
     &       ' IT MUST BE GREATER THAN 0.'//)

99999 CONTINUE

#ifdef CMPI
      call msg_fini()
#endif
      WRITE(*,'(a)') 'RADIATION STRESS READ ERROR (RS2-1)'
      CALL EXIT(1)

      END SUBROUTINE


!***********************************************************************
!   SUBROUTINE RS2GET
!***********************************************************************

      subroutine rs2get(rsnx,rsny,np)
      implicit none

      integer,intent(in) :: np
      real(kind=8),intent(out) :: rsnx(:),rsny(:)

      if(type.eq.1) then
        write(*,'(a)') 'RS2GET: SORRY! TYPE-1 RS FILE NOT YET SUPPORTED.'
        CALL EXIT(1) 
!        call rs2get_type1(rsnx,rsny,np)
      else if(type.eq.2) then
        call rs2get_type2(rsnx,rsny,np)
      else
         write(*,'(a)') 'RS2 TYPE IS NEITHER 1 NOR 2.  TERMINATED.'
         CALL EXIT(1) 
      endif

      return

      end subroutine

!$$$C***********************************************************************
!$$$C   SUBROUTINE RS2GET_TYPE1
!$$$C***********************************************************************
!$$$
!$$$      subroutine rs2get_type1(rsnx,rsny,np)
!$$$      use sizes
!$$$#ifdef CMPI     
!$$$      use messenger, only: msg_fini, msg_ibcast, msg_cbcast, msg_rbcast
!$$$#ifdef HAVE_MPI_MOD
!$$$      use mpi
!$$$      implicit none
!$$$#else
!$$$      implicit none
!$$$      include 'mpif.h'
!$$$#endif
!$$$#else
!$$$      implicit none
!$$$#endif 
!$$$      integer,intent(in) :: np
!$$$      real(kind=8),intent(out) :: rsnx(:),rsny(:)
!$$$      integer :: i,j,k,s,n,ni,nj,fid,ios(1)
!$$$      integer :: isActive(1)
!$$$      character(100) dummyc, spectrum_id(1)
!$$$      real(kind=8), dimension(:,:,:), allocatable :: rsxy
!$$$      real(kind=8) :: rsx,rsy,rsx1,rsx2,rsx3,rsy1,rsy2,rsy3,cnf1,cnf2,cnf3
!$$$      real(kind=8) :: angle_st,pi
!$$$
!$$$      write(*,'(a)') 'RS2GET_TYPE1: SORRY! TYPE-1 RS FILE NOT YET SUPPORTED.'
!$$$      CALL EXIT(1)
!$$$
!$$$      if(myproc.eq.0) then
!$$$        write(screenunit,*)
!$$$        write(screenunit,*) 'READING IN RADIATION STRESS'
!$$$      endif
!$$$      write(16,*)
!$$$      write(16,*) 'READING IN RADIATION STRESS'
!$$$
!$$$      pi = 4.e0*atan(1.e0)
!$$$
!$$$      !increment counter
!$$$      cntSnaps = cntSnaps+1
!$$$
!$$$      !initialize
!$$$      do i=1,np
!$$$        rsnx(I)=0.e0
!$$$        rsny(I)=0.e0
!$$$      enddo
!$$$
!$$$
!$$$      fid = 23
!$$$#ifdef CMPI
!$$$      if(myproc.eq.0) read(fid,'(A)',iostat=ios(1)) dummyc
!$$$      call msg_ibcast(ios,1)
!$$$      call msg_cbcast(dummyc,100)
!$$$#else
!$$$      read(fid,'(A)',iostat=ios(1)) dummyc
!$$$#endif
!$$$      write(16,*) 'READING RADIATION STRESS SNAP.  '//
!$$$     &              'SNAP ID = ', trim(dummyc)
!$$$      !if fort.23 had reached EOF or got any read error,
!$$$      ! stop reading fort.23 and use insert blank snaps from now on.
!$$$      if(ios(1) < 0) then
!$$$        if(myproc.eq.0) then
!$$$          write(screenunit,*)
!$$$          write(screenunit,*) 'UNIT23 FILE RUN OUT.'
!$$$          write(screenunit,*) 'BLANK SNAP INSERTED FROM NOW ON.'
!$$$        endif
!$$$        write(16,*)
!$$$        write(16,*) 'UNIT23 FILE RUN OUT.'
!$$$        write(16,*) 'BLANK SNAP WILL BE INSERTED FROM NOW ON.'
!$$$        do s=1,ns
!$$$          endoffile(s) = 1
!$$$        enddo
!$$$        return
!$$$      endif
!$$$
!$$$!-----------
!$$$! S-Loop
!$$$!-----------
!$$$
!$$$      do s=1,ns   ! put a blank snap for the first nbs(s) snaps
!$$$        if(cntsnaps.le.nbs(s)) then
!$$$          if(abs(nscreen) >= 1) then
!$$$            if(myproc.eq.0) then
!$$$              write(screenunit,*) 'INSERTING A BLANK RS SNAP, REGION = ',s,
!$$$     &                   ', COUNT = ',cntsnaps
!$$$            endif
!$$$            write(16,*) 'INSERTING A BLANK RS SNAP, REGION = ',s,
!$$$     &                  ', COUNT = ',cntsnaps
!$$$          endif
!$$$          cycle    
!$$$        endif
!$$$
!$$$        ni = nis(s)
!$$$        nj = njs(s)
!$$$
!$$$        allocate(rsxy(2,nis(s),njs(s)))
!$$$
!$$$        fid = 23
!$$$        if (myproc.eq.0) read(fid,*,iostat=ios(1)) ! skip a line
!$$$#ifdef CMPI     
!$$$        call msg_ibcast(ios,1)
!$$$#endif
!$$$        if (ios(1) > 0) g.oto 9999
!$$$
!$$$        if(myproc.eq.0) read(fid,*,iostat=ios(1)) isActive(1) ! active flag
!$$$
!$$$#ifdef CMPI        
!$$$        call msg_ibcast(ios,1)
!$$$        call msg_ibcast(isActive,1)
!$$$#endif
!$$$
!$$$        if(ios(1) > 0)       g.oto 9999
!$$$        if(isActive(1).eq.0) g.oto 100
!$$$
!$$$        if(myproc.eq.0) then
!$$$          do j = nj, 1, -1
!$$$            read (fid,'(5E15.7)',iostat=ios(1))
!$$$     &           ((rsxy(k,i,j),k=1,2),i=1,ni)
!$$$          enddo
!$$$        endif
!$$$
!$$$#ifdef CMPI      
!$$$        call msg_ibcast(ios,1)
!$$$        call msg_rbcast(rsxy,2*ni*nj)
!$$$#endif
!$$$        if(ios(1) > 0)       g.oto 9999
!$$$
!$$$C--------------------------------------------------------------------------
!$$$C   Interpolate RS Snap onto local ADCIRC mesh
!$$$C--------------------------------------------------------------------------
!$$$
!$$$        do n=1,np
!$$$          if(rsij(1,n).ne.s) cycle !node n does not reside in STWAVE grid s
!$$$            
!$$$          i = rsij(3,n)
!$$$          j = rsij(4,n)
!$$$
!$$$          if(rsij(2,n).eq.1) then !node n resides in triangle 1 in cell (i,j).
!$$$            rsx1 = rsxy(1,i,j)
!$$$            rsy1 = rsxy(2,i,j)
!$$$            rsx2 = rsxy(1,i+1,j)
!$$$            rsy2 = rsxy(2,i+1,j)
!$$$            rsx3 = rsxy(1,i+1,j+1)
!$$$            rsy3 = rsxy(2,i+1,j+1)
!$$$          else                    !node n resides in triangle 2 in cell (i,j).
!$$$            rsx1 = rsxy(1,i+1,j+1)
!$$$            rsy1 = rsxy(2,i+1,j+1)
!$$$            rsx2 = rsxy(1,i,j+1)
!$$$            rsy2 = rsxy(2,i,j+1)
!$$$            rsx3 = rsxy(1,i,j)
!$$$            rsy3 = rsxy(2,i,j)
!$$$          endif
!$$$
!$$$          !interpolation
!$$$          rsx= rscnf(1,n)*rsx1 + rscnf(2,n)*rsx2 + rscnf(3,n)*rsx3
!$$$          rsy= rscnf(1,n)*rsy1 + rscnf(2,n)*rsy2 + rscnf(3,n)*rsy3
!$$$
!$$$          !need to rotate
!$$$          angle_st = (azimuths(s)*pi)/180.e0
!$$$          rsnx(n) = cos(angle_st)*rsx-sin(angle_st)*rsy
!$$$          rsny(n) = sin(angle_st)*rsx+cos(angle_st)*rsy
!$$$        enddo
!$$$
!$$$ 100    continue
!$$$        deallocate(rsxy)
!$$$
!$$$      enddo       ! end s-loop
!$$$
!$$$      
!$$$      if(myproc.eq.0) then
!$$$        write(screenunit,*) 'RADIATION STRESS IS LOADED'
!$$$        write(screenunit,*)
!$$$      endif
!$$$      write(16,*) 'RADIATION STRESS IS LOADED'
!$$$      write(16,*)
!$$$
!$$$      return
!$$$
!$$$ 9999 continue
!$$$      if(myproc.eq.0) then
!$$$        write(screenunit,*) trim(radfiles(s)),' IS NOT COMPLETE.'
!$$$        write(screenunit,*) 'PROGRAM WILL BE TERMINATED.'
!$$$      endif
!$$$      write(16,*) trim(radfiles(s)),' IS NOT COMPLETE.'
!$$$      write(16,*) 'PROGRAM WILL BE TERMINATED.'
!$$$      close(16)
!$$$
!$$$#ifdef CMPI
!$$$      call msg_fini()
!$$$#endif
!$$$      CALL EXIT(1)
!$$$
!$$$      end subroutine

!***********************************************************************
!   SUBROUTINE RS2GET_TYPE2
!***********************************************************************

      subroutine rs2get_type2(rsnx,rsny,np)
      use sizes
#ifdef CMPI
!asey 090327: Implement Seizo's changes for buffering the radiation stress gradients.     
      use messenger, only: msg_fini, msg_ibcast, msg_cbcast, msg_rbcastd
#ifdef HAVE_MPI_MOD
#ifndef ROSECOMP
      use mpi
#endif
      implicit none
#else
      implicit none
      include 'mpif.h'
#endif 
#else
      implicit none
#endif 

      integer,intent(in) :: np
      real(kind=8),intent(out) :: rsnx(:),rsny(:)
      integer :: i,j,k,s,n,ni,nj,fid,ios(1)
      integer :: isActive(1)
      character(100) dummyc, spectrum_id(1)
      real(kind=8), dimension(:,:,:), allocatable :: rsxy
      real(kind=8) :: rsx,rsy,rsx1,rsx2,rsx3,rsy1,rsy2,rsy3,cnf1,cnf2,cnf3
      real(kind=8) :: angle_st,pi

      if(myproc.eq.0) then
        write(screenunit,*)
        write(screenunit,*) 'READING IN RADIATION STRESS'
      endif
      write(16,*)
      write(16,*) 'READING IN RADIATION STRESS'

      pi = 4.e0*atan(1.e0)

      !increment counter
      cntSnaps = cntSnaps+1

      !initialize
      do i=1,np
        rsnx(I)=0.e0
        rsny(I)=0.e0
      enddo


!-----------
! S-Loop
!-----------

      do s=1,ns   ! put a blank snap for the first nbs(s) snaps
        if(cntsnaps.le.nbs(s)) then
          if(abs(nscreen) >= 1) then
            if(myproc.eq.0) then
         write(screenunit,*) 'INSERTING A BLANK RS SNAP, REGION = ',s,
     &                   ', COUNT = ',cntsnaps
            endif
            write(16,*) 'INSERTING A BLANK RS SNAP, REGION = ',s,
     &                  ', COUNT = ',cntsnaps
          endif
          cycle    
        endif

        ! if previous EOF write blank RS SNAP
        if(endoffile(s).eq.1) then
          if (myproc.eq.0) then
             write(screenunit,*) trim(radfiles(s)),
     &                          ' PREVIOUSLY REACHED EOF.'
             write(screenunit,'(A25,I6,A9,I6)') 
     &         " BLANK RS SNAP, REGION = ",s, " COUNT = ", cntsnaps
          endif
          cycle
        endif          

        ni = nis(s)
        nj = njs(s)

        fid = rfids(s)

        if (endoffile(s) .eq. 0) then
          if (myproc.eq.0) then
            read (fid,'(A)',iostat=ios(1)) spectrum_id(1) ! spectrum identifier
          endif
#ifdef CMPI    
          call msg_ibcast(ios,1)
#endif

          if (ios(1) < 0) then
            endoffile(s) = 1
            if (myproc.eq.0) then
              write(screenunit,*) trim(radfiles(s)),' REACHED EOF.'
              write(screenunit,'(A25,I6,A9,I6)') 
     &           " BLANK RS SNAP, REGION = ",s, " COUNT = ", cntsnaps
            endif
            cycle
          else
            if (myproc.eq.0) then
        write(screenunit,'(A)') " spectrum_id = "//trim(spectrum_id(1))
              write(screenunit,'(A19,I6,A9,I6)') 
     &         " RS SNAP, REGION = ",s, " COUNT = ", cntsnaps
            endif
          endif
        endif
         
        allocate(rsxy(2,nis(s),njs(s)))

#ifdef CMPI     
        if(myproc.eq.0) then
          ios(1) = 0
          do while(.true.)
            do j = nj, 1, -1
              read (fid,*,err=511) 
     &           ((rsxy(k,i,j),k=1,2),i=1,ni)
            enddo
            exit
511         ios(1) = 1
            exit
          enddo
        endif
        call msg_ibcast(ios,1)
        if(ios(1) == 0) then        
!asey 090327: Implement Seizo's changes for buffering the radiation stress gradients.
          call msg_rbcastd(rsxy,2,ni,nj)
        else
          call rs2get_type2_terminate(radfiles(s))
        endif
#else
        do j = nj, 1, -1
          read (fid,*,err=9999) 
     &         ((rsxy(k,i,j),k=1,2),i=1,ni)
        enddo
#endif

!--------------------------------------------------------------------------
!   Interpolate RS Snap onto local ADCIRC mesh
!--------------------------------------------------------------------------

        do n=1,np
          if(rsij(1,n,s).eq.0) cycle !node n does not reside in STWAVE grid s
            
          i = rsij(2,n,s)
          j = rsij(3,n,s)

          if(rsij(1,n,s).eq.1) then !node n resides in triangle 1 in cell (i,j).
            rsx1 = rsxy(1,i,j)
            rsy1 = rsxy(2,i,j)
            rsx2 = rsxy(1,i+1,j)
            rsy2 = rsxy(2,i+1,j)
            rsx3 = rsxy(1,i+1,j+1)
            rsy3 = rsxy(2,i+1,j+1)
          else                    !node n resides in triangle 2 in cell (i,j).
            rsx1 = rsxy(1,i+1,j+1)
            rsy1 = rsxy(2,i+1,j+1)
            rsx2 = rsxy(1,i,j+1)
            rsy2 = rsxy(2,i,j+1)
            rsx3 = rsxy(1,i,j)
            rsy3 = rsxy(2,i,j)
          endif

          !interpolation
        rsx= rscnf(1,n,s)*rsx1 + rscnf(2,n,s)*rsx2 + rscnf(3,n,s)*rsx3
        rsy= rscnf(1,n,s)*rsy1 + rscnf(2,n,s)*rsy2 + rscnf(3,n,s)*rsy3

       if((rsnx(n)*rsnx(n)+rsny(n)*rsny(n)).lt.(rsx*rsx+rsy*rsy))then
            !need to rotate
            angle_st = (azimuths(s)*pi)/180.e0
            rsnx(n) = cos(angle_st)*rsx-sin(angle_st)*rsy
            rsny(n) = sin(angle_st)*rsx+cos(angle_st)*rsy
          endif
        enddo

        deallocate(rsxy)

      enddo       ! end s-loop

      if(myproc.eq.0) then
        write(screenunit,*) 'RADIATION STRESS IS LOADED'
        write(screenunit,*)
      endif
      write(16,*) 'RADIATION STRESS IS LOADED'
      write(16,*)

      return

 9999 call rs2get_type2_terminate(radfiles(s))

      end subroutine rs2get_type2

      subroutine rs2get_type2_terminate(filename)
        implicit none
        character(*),intent(in) :: filename
        if(myproc.eq.0) then
          write(screenunit,*) trim(filename),' IS NOT COMPLETE.'
          write(screenunit,*) 'PROGRAM WILL BE TERMINATED.'
        endif
        write(16,*) trim(filename),' IS NOT COMPLETE.'
        write(16,*) 'PROGRAM WILL BE TERMINATED.'
        close(16)
#ifdef CMPI
        call msg_fini()
#endif
        CALL EXIT(1)
      end subroutine rs2get_type2_terminate

    
      END MODULE RS2
      
      
