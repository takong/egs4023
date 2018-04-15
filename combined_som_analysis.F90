        Program SOM
        
        include 'netcdf.inc'
        parameter(ni=100,nj=100, ndata=1000)
        character*40 InFile
        character*13 vname
        real lon(ni),lat(nj)
        real lon2D(ni,nj), lat2D(ni,nj)
        character*40 OutFile1, OutFile
        character*13 tname !vname(nvar),
        character*40 TimeUnitString
        character*100 command
        character*60 datafile
        character*5 TheLon, TheLat
        character*40 InFile(nvar)
        character*100 ref_date,time_calendar
        real time(ntime)
        real ExtrDays(nwarea,ntime)
        dimension ilev(nvar)
        real var(ni,nj)

        !------Open Input file
        write(*,*) InFile(1)

        iret = nf_open(InFile(1),NF_NOWRITE,ncid1)
        if(iret .ne. NF_NOERR) call handle_error(iret)
        call ReadLonLatTime(lon,lat,time,ni,nj,nk,itime_start,ntime,ncid1)

        !------Close the file
        iret = nf_close(ncid1)
        if (ret.ne.NF_NOERR) then
        call handle_error (iret)
        endif

        !------------------------------------------------------------
        !     Set the study domain (i.e. Southern Africa)
        !------------------------------------------------------------
        call Set_Domain_Ids(lon, lat, ni, nj, &
        rlonD_start, rlonD_stop,rlatD_start, rlatD_stop, &
        lonD_start,  lonD_stop, latD_start,  latD_stop )
        nlonD = lonD_stop - lonD_start +1
        nlatD = latD_stop - latD_start +1
        !------------------------------------------------
        !    Prepare and open the output file (NC)
        !-------------------------------------------------
        !     if (idata .eq. 1) call Prepare_CombinedSOMsOutputFile(OutFile1, nlonD, nlatD, &
        !                   lonD_start,  lonD_stop, latD_start,  latD_stop, &
        !                   lon, lat, ni, nj,time, ntime, ndata, ref_date, time_calendar)

        call Prepare_OutputFile(OutFile, nlonD, nlatD, &
        lonD_start,  lonD_stop, latD_start,  latD_stop, &
        lon, lat, ni, nj,time, ntime, nclass,nwarea, ref_date, time_calendar)

        write(*,*)"Done creating ", trim(OutFile)


        !----------------------------------------------
        !     Read and write the neccessary variables
        !----------------------------------------------
        do ivar = 1, nuse
        write(*,*) vname(ivar)
        Call ReadWriteAllVars(OutFile, InFile(ivar), ilev(ivar), vname(ivar), nlonD, nlatD, &
        lonD_start, latD_start, itime_start, ntime)
        write(*,*) "done for ", OutFile, ":", vname(ivar), nuse
        enddo

        !---------------------------------
        !      Perform SOMs analysis
        !---------------------------------
        call GenSOMs(OutFile, InFile(1),ilev(1), vname(1), nlonD, nlatD, lonD_start, &   
        latD_start, ExtrDays, itime_start, ntime, idata, nwarea)
        
        write(*,*) "Analysis successfully completed !!!"

        end program SOM

!======================================================================
!======================================================================
      Subroutine GenSOMs(OutFile,InFile,ilev,vname, ni,nj,lon1,lat1,ExtrDays,itime_start,ntime,idata, nwarea)
      include 'netcdf.inc'
      character*100 command
      character*40 OutFile
      character*40 InFile
      character*13 vname
      real var(ni,nj),ExtrDays(nwarea,ntime)
      character*15 TheDate(ntime)
      dimension pcvar(ni*nj), iDays(ntime)


        open(1, file= "x.asc")
        idt = 0
        nprow = 0

        open(8, file="iused_idx.dat")

!------Open Input file and read data
          iret = nf_open(InFile,NF_NOWRITE,ncid1)
          if(iret .ne. NF_NOERR) call handle_error(iret)
        iwarea = 6
        do  itime = 1,ntime           
         if (ExtrDays(iwarea,itime) .eq. 1) then
           idt = idt + 1
           itread = itime + (itime_start-1)
            call Read_RAWData(vname,var,ni, nj, ilev, lon1, lat1, itread, ncid1)
            if (trim(vname) .eq. 'psl') var = var/100.
            if (trim(vname) .eq. 'PRM') then
            var =  var*86400.0
             do i = 1, ni
              do j = 1, nj
               if (var(i,j) .gt. 1000.00) var(i,j) = 0.0
              enddo
             enddo
            endif

            npcol = 0   
            do j= 1,nj
              do i = 1,ni
               npcol = npcol + 1
               pcvar(npcol) = var(i,j)
             enddo
            enddo
!            write(*,*) "idt:",  idt, npcol
            if(idt .eq. 1) then        
              write(*,*) "npcol :", npcol  
              write(1,*) npcol 
              if (idata .eq. 1)  write(11,*) npcol
              npcol1 = npcol
            endif
            if (npcol .eq. npcol1) then
              write(1, 10) (pcvar(ip),ip=1,npcol)
              write(11,10) (pcvar(ip),ip=1,npcol)
              nprow = nprow + 1 
              write(8,*) itime
              write(88,*) idata, itime
            else
              write(*,*)"inconsitent with no of column...", idt
            endif 
          endif 
         enddo
         write(*,*)"Total number of days used:", nprow, ntime
         write(22,*) idata, npcol, nprow
        close(8)
        close(1)

!===== close file
        iret = nf_close(ncid1)
        if (ret.ne.NF_NOERR) then
        call handle_error (iret)
        endif

  10   format(100000(f8.2))
  11   format((f10.2))

       
!---------------------------------------------------------------
!(3)  Run the SOMS Codes
!---------------------------------------------------------------
      nxdim = 3
      nydim = 4
      if(nprow > 0) then
!----Initialize the data saying the som will be 4X4 nodes using a rectablular topology
    write(command,*)"./randinit -din x.asc -cout x.cod -xdim", nxdim, " -ydim ", nydim, " -topol rect -neigh bubble"
       write(*,*) command
       call system(command)

!----This does the SOM in 2 phases, one a course pass with 50000 iterations and if 
!    you want a 2nd one use a much higher number of iterations than the first
       write(command,*)"./vsom -din x.asc -cin x.cod -cout xI.cod -rlen 50000 -alpha 0.1 -radius 2"
       write(*,*) command
       call system(command)

       write(command,*)"./vsom -din x.asc -cin xI.cod -cout xII.cod -rlen 50000 -alpha 0.01 -radius 1"
       write(*,*) command
       call system(command)

!---This produces the Sammon maps to show how the nodes relate to each other in data space
       write(command,*)"./sammon -cin xI.cod -cout xI.sam -rlen 1000 -ps"
       write(*,*) command
       call system(command)

       write(command,*)"./sammon -cin xII.cod -cout xII.sam -rlen 1000 -ps"
       write(*,*) command
       call system(command)

!   Makes Frequency Map
       write(command,*)"./visual -din  x.asc -cin xII.cod -dout xII.vis"
       write(*,*) command
       call system(command)

       write(command,*)"cp  xII.vis  ", trim(OutFile) // ".vis"
       write(*,*) command
       call system(command)

       write(cammand,*)"./annfreq xII"
       write(*,*) command
       call system(command)

!----Get Archetype Matrix, run on linux
       write(command,*)"./cod2grads xII"
       write(*,*) command
       call system(command)
      else
          write(*,*)'No pca run because number of rows is:', nprow 
      endif


!---------------------------------------------------------------
!(4)  Read & Write SOMs Data
!---------------------------------------------------------------
!------ Open output file
       iret = nf_open(OutFile,NF_WRITE,ncido)
       if(iret .ne. NF_NOERR) call handle_error(iret)

       npfac = nxdim * nydim
       Call Read_SOMSData(npcol,npfac,nxdim,nydim,ni,nj,ncido)

!---------------------------------------------------------------
! Close output file
!--------------------------------------------------------------
       iret = nf_close(ncido)
       if (ret.ne.NF_NOERR) then
       call handle_error (iret)
       endif

      return
      end

      Subroutine GenCombinedSOMs(OutFile, ni,nj,ndata)
      include 'netcdf.inc'
      character*100 command
      character*40 OutFile
      character*40 InFile
      character*13 vname
      real var(ni,nj)
      dimension pcvar(ni*nj)

      call system("cp all_soms_x.asc ax.asc")

      open(11, file= 'ax.asc')
      read(11,*) npcol
      write(*,*) npcol
      close(11)


!       go to 999

      nprow =1
!---------------------------------------------------------------
!(3)  Run the SOMS Codes
!---------------------------------------------------------------
      nxdim = 4
      nydim = 3
      if(nprow > 0) then
!----Initialize the data saying the som will be 4X4 nodes using a rectablular topology
       write(command,*)"./randinit -din x.asc -cout x.cod -xdim", nxdim, " -ydim ", nydim, " -topol rect -neigh bubble"
       write(*,*) command
       call system(command)

!----This does the SOM in 2 phases, one a course pass with 50000 iterations and if 
!    you want a 2nd one use a much higher number of iterations than the first
       write(command,*)"./vsom -din ax.asc -cin x.cod -cout xI.cod -rlen 50000 -alpha 0.1 -radius 3"
       write(*,*) command
       call system(command)

       write(command,*)"./vsom -din ax.asc -cin xI.cod -cout xII.cod -rlen 645000 -alpha 0.01 -radius 1"
       write(*,*) command
       call system(command)

!---This produces the Sammon maps to show how the nodes relate to each other in data space
       write(command,*)"./sammon -cin xI.cod -cout xICom.sam -rlen 1000 -ps"
       write(*,*) command
       call system(command)

       write(command,*)"./sammon -cin xII.cod -cout xIICom.sam -rlen 1000 -ps"
       write(*,*) command
       call system(command)

!   Makes Frequency Map
       write(command,*)"./visual -din  ax.asc -cin xII.cod -dout xII.vis"
       write(*,*) command
       call system(command)

       write(command,*)"cp  xII.vis  ", trim(OutFile) // "allsoms.vis"
       write(*,*) command
       call system(command)

       write(cammand,*)"./annfreq xII"
       write(*,*) command
       call system(command)

!----Get Archetype Matrix, run on linux
       write(command,*)"./cod2grads xII"
       write(*,*) command
       call system(command)
      else
          write(*,*)'No pca run because number of rows is:', nprow 
      endif

      999 continue
!---------------------------------------------------------------
!(4)  Read & Write SOMs Data
!---------------------------------------------------------------
!------ Open output file
       iret = nf_open(OutFile,NF_WRITE,ncido)
       if(iret .ne. NF_NOERR) call handle_error(iret)
       write(*,*) OutFile
       npfac = nxdim * nydim
       Call Read_CombinedSOMSData(npcol,npfac,nxdim,nydim,ni,nj,ncido,ndata)

!---------------------------------------------------------------
! Close output file
!--------------------------------------------------------------
       iret = nf_close(ncido)
       if (iret.ne.NF_NOERR) then
       call handle_error (iret)
       endif
       write(*,*)"bye!"

      return
      end
!=====================================================================
       Subroutine Read_CombinedSOMSData(npcol,npfac,nxdim,nydim,nlon,nlat,ncido,ndata)
 
       real pcvar(npcol,npfac),Freq(npfac)
       real pcvar2d(nlon,nlat), SubFreq(ndata,npfac), Freq1D(ndata)
       real DataAmount(ndata),  DataStop(ndata)
       character*13 vname
       
       write(*,*) npcol,npfac,nlon,nlat
!===== Read results PCA output (loadings)      
       open(1, file = "xII.cod")   
        read (1,*)  ! npcol
          do ifac = 1,npfac
           read(1,*) (pcvar(icol,  ifac),icol=1,npcol)
!           write(*,*) (pcvar(icol,ifac),icol=1,npcol)
          enddo 
       close(1)


      close(88)

       open(11, file="all_soms_info")

       LastData = 0

       write(*,*) "ndata", ndata
       do idata = 1,ndata
         read(11,*) tmp1, tmp2, DataAmount(idata)
         DataStop(idata) = LastData +  DataAmount(idata)
         LastData = DataStop(idata)
       enddo


!===== Calculate the Frequecy
       Freq(:) = 0
       tFreq   = 0
       icount   = 0 
       SubFreq(:,:) = 0.0
       open(1, file = "xII.vis")
       open(88, file="all_iused_idx.dat")
       read (1, *,end=999)
        do while (.true.)
         read (1, *, end=999) ifxin,ifyin
         write(*,*) ifxin,ifyin
         ifacin = ifyin*nxdim + ifxin + 1
         Freq(ifacin) = Freq(ifacin) + 1   
         tFreq = tFreq + 1 
         icount = icount + 1
         idin = 1

!         do idata = 2,ndata
!          if (DataStop(idata-1) < icount) idin = idata   
!         enddo

         read(88,*) idata, iused 
         write(*,*) idata, iused,ifacin 
         SubFreq(idata,ifacin) = SubFreq(idata,ifacin) + 1
         call Write_CombinedSOMDataTime(ifacin,idata,iused,ncido)       
       enddo

  999  close(1)
       close(88)


        do ifac = 1,npfac
         fsum = 0
         do idata = 1,ndata
           fsum = fsum + SubFreq(idata,ifac)
         enddo
         write(*,*) fsum
         do idata = 1,ndata
            SubFreq(idata,ifac) = SubFreq(idata,ifac) !/fsum * 100.0
         enddo
        enddo

       do ifac = 1,npfac
!===== Convert SOMs Output (nodes) to 2D data
         icol = 0
         do j = 1,nlat
          do i = 1,nlon
           pcvar2d(i,j) = -1.0e+30
           icol = icol + 1
             pcvar2d(i,j) = pcvar(icol,ifac)
          enddo
         enddo

!====== Output SOMs Data
        write(*,*)ifac,Freq(ifac),tFreq 
        iFrq = Freq(ifac)/tFreq*10000 
        Frq = iFrq/100.
        Freq1D = SubFreq(:,ifac)
        call Write_CombinedSOMData(pcvar2d,Frq,Freq1D,nlon,nlat,ndata,ifac,ncido)
      enddo
      end

!======================================================================
      Subroutine Write_CombinedSOMData(var,Freq,Freq1D,ni,nj,ndata,ipfac,ncido)
      include 'netcdf.inc'
      real var(ni,nj), Freq1D(ndata)
      character*13 vname
      integer data_start(3),data_count(3)
      integer data_start2(2),data_count2(2)


!        write(*,*)"1", ni,nj,Freq,ipfac,ndata
 
        data_start(1) = 1
        data_start(2) = 1
        data_start(3) = ipfac
                
        data_count(1) = ni
        data_count(2) = nj
        data_count(3) = 1 
        iret = nf_inq_varid (ncido, "SOMS_Nodes", ivarid)
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_put_vara_real (ncido,ivarid,data_start,data_count,var)
        if (iret.ne.NF_NOERR) call handle_error(iret)

!        write(*,*)"2", ni,nj,Freq,ipfac,ndata
        data_start2(1) = 1
        data_start2(2) = ipfac
                
        data_count2(1) = ndata
        data_count2(2) = 1 
        iret = nf_inq_varid (ncido, "SubFreq", ivarid)
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_put_vara_real (ncido,ivarid,data_start2,data_count2,Freq1D)
        if (iret.ne.NF_NOERR) call handle_error(iret)

!        write(*,*)"3", ni,nj,Freq,ipfac,ndata
        vfac =ipfac
        idata_start = ipfac
        idata_count = 1 
        iret = nf_inq_varid (ncido, 'factor', ivarid)
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_put_vara_real (ncido,ivarid,idata_start,idata_count,vfac)
        if (iret.ne.NF_NOERR) call handle_error(iret)

!        write(*,*)"4", ni,nj,Freq,ipfac,ndata
        vfac = Freq
        idata_start = ipfac
        idata_count = 1 
        iret = nf_inq_varid (ncido, 'AllFreq', ivarid)
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_put_vara_real (ncido,ivarid,idata_start,idata_count,vfac)
        if (iret.ne.NF_NOERR) call handle_error(iret)


      return
      end


!=========================================================================
      Subroutine Write_CombinedSOMDataTime(ipfac,idata,itime,ncido)
      include 'netcdf.inc'
      integer idata_start(2),idata_count(2)

        write(*,*) idata,itime,ipfac
        vfac = ipfac

        idata_start(1) = idata
        idata_start(2) = itime

        idata_count(1) = 1 
        idata_count(2) = 1 

        iret = nf_inq_varid (ncido, 'NodeTime', ivarid)
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_put_vara_real (ncido,ivarid,idata_start,idata_count,vfac)
        if (iret.ne.NF_NOERR) call handle_error(iret)

      return
      end

!=====================================================================
       Subroutine Read_SOMSData(npcol,npfac,nxdim,nydim,nlon,nlat,ncido)
       real pcvar(npcol,npfac),Freq(npfac)
       real pcvar2d(nlon,nlat)
       character*13 vname

       write(*,*) npcol,npfac,nlon,nlat

!===== Read results PCA output (loadings)      
       open(1, file = "xII.cod")   
        read (1,*)  ! npcol
          do ifac = 1,npfac
           read(1,*) (pcvar(icol,  ifac),icol=1,npcol)
!           write(*,*) (pcvar(icol,ifac),icol=1,npcol)
          enddo 
       close(1)

!===== Calculate the Frequecy
       Freq(:) = 0
       tFreq  = 0
       itime = 0
       open(1, file = "xII.vis")
       open(8, file="iused_idx.dat")
       read (1, *,end=999)
        do while (.true.)
         itime = itime + 1
         read (1, *, end=999) ifxin,ifyin
         ifacin = ifyin*nxdim + ifxin + 1
         Freq(ifacin) = Freq(ifacin) + 1   
         tFreq = tFreq + 1                 
         if(ifacin .gt. npfac) stop 'ifac_error'
         read(8,*) iused 
         call Write_SOMDataTime(ifacin,iused,ncido)
       enddo
  999  close(1)
       close(8)


       do ifac = 1,npfac
!===== Convert SOMs Output (nodes) to 2D data
         icol = 0
         do j = 1,nlat
          do i = 1,nlon
           pcvar2d(i,j) = -1.0e+30
           icol = icol + 1
             pcvar2d(i,j) = pcvar(icol,ifac)
          enddo
         enddo

!====== Output SOMs Data
        write(*,*)ifac,Freq(ifac),tFreq 
        iFrq = Freq(ifac)/tFreq*10000 
        Frq = iFrq/100.
        call Write_SOMData(pcvar2d,Frq,nlon,nlat,ifac,ncido)
      enddo

      end



!======================================================================
      Subroutine Write_SOMData(var,Freq,ni,nj,ipfac,ncido)
      include 'netcdf.inc'
      real var(ni,nj)
      character*13 vname
      integer data_start(3),data_count(3)

        data_start(1) = 1
        data_start(2) = 1
        data_start(3) = ipfac
       
        data_count(1) = ni
        data_count(2) = nj
        data_count(3) = 1 
        iret = nf_inq_varid (ncido, "SOMS_Nodes", ivarid)
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_put_vara_real (ncido,ivarid,data_start,data_count,var)
        if (iret.ne.NF_NOERR) call handle_error(iret)

        vfac =ipfac
        idata_start = ipfac
        idata_count = 1 
        iret = nf_inq_varid (ncido, 'factor', ivarid)
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_put_vara_real (ncido,ivarid,idata_start,idata_count,vfac)
        if (iret.ne.NF_NOERR) call handle_error(iret)

        vfac = Freq
        idata_start = ipfac
        idata_count = 1 
        iret = nf_inq_varid (ncido, 'Freq', ivarid)
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_put_vara_real (ncido,ivarid,idata_start,idata_count,vfac)
        if (iret.ne.NF_NOERR) call handle_error(iret)


      return
      end



!=========================================================================
      Subroutine Write_SOMDataTime(ipfac,itime,ncido)
      include 'netcdf.inc'

        write(*,*) itime, ipfac
        vfac =ipfac
        idata_start = itime
        idata_count = 1 
        iret = nf_inq_varid (ncido, 'NodeTime', ivarid)
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_put_vara_real (ncido,ivarid,idata_start,idata_count,vfac)
        if (iret.ne.NF_NOERR) call handle_error(iret)

      return
      end



!=========================================================================
      Subroutine Read_RAWData(vname,var,nlon, nlat, ilev, lon1, lat1, itime, ncid1)
                              
      include 'netcdf.inc'
      real var(nlon,nlat)
      character*13 vname
      integer data3_start(3),data3_count(3)
      integer data4_start(4),data4_count(4)

        if (ilev .eq. 0) then
          data3_start(1) = lon1
          data3_start(2) = lat1
          data3_start(3) = itime

          data3_count(1) = nlon
          data3_count(2) = nlat
          data3_count(3) = 1 
          iret = nf_inq_varid (ncid1, vname, ivarid)
          if (iret.ne.NF_NOERR) call handle_error (iret)
          iret = nf_get_vara_real (ncid1,ivarid,data3_start,data3_count,var)
        else
          data4_start(1) = lon1
          data4_start(2) = lat1
          data4_start(3) = ilev
          data4_start(4) = itime

          data4_count(1) = nlon
          data4_count(2) = nlat
          data4_count(3) = 1
          data4_count(4) = 1 
          iret = nf_inq_varid (ncid1, vname, ivarid)
          if (iret.ne.NF_NOERR) call handle_error (iret)
          iret = nf_get_vara_real (ncid1,ivarid,data4_start,data4_count,var)
        endif
      return
      end


!================================================================
      Subroutine Set_Domain_Ids(lon, lat, ni, nj, &
          rlon_start, rlon_stop, rlat_start, rlat_stop, &
           lon_start,  lon_stop,  lat_start,  lat_stop )

      real lon(ni),lat(nj)

!----------- Set lon & lat indexes (i,j) to search region of interest
         lon_start = 1
         lon_stop  = ni
         lat_start = 1
         lat_stop  = nj

        do i=1,ni
         if(lon(i)<=rlon_start .and. lon(i+1)>=rlon_start) lon_start=i
         if(lon(i)<=rlon_stop  .and. lon(i+1)>=rlon_stop)  lon_stop=i
        enddo  
        do j=1,nj
         if(lat(j)<=rlat_start .and. lat(j+1)>=rlat_start) lat_start=j
         if(lat(j)<=rlat_stop  .and. lat(j+1)>=rlat_stop ) lat_stop=j

!         if(lat(j+1)<=dlat_start .and. lat(j)>=dlat_start) lat_start=j
!         if(lat(j+1)<=dlat_stop  .and. lat(j)>=dlat_stop ) lat_stop=j

        enddo
         if(lat_start > lat_stop) then
            mtmp = lat_start
            lat_start = lat_stop
            lat_stop = mtmp
         endif

      return
      end

!==================================================================
      Subroutine Prepare_OutputFile(OutFile, nlon, nlat, &
                   lon_start,  lon_stop, lat_start,  lat_stop, &
                   lon, lat, ni, nj,time, ntime,nclass, nwarea, ref_date, time_calendar)

      include 'netcdf.inc'
      character*40 OutFile
      character*100 command
      real lon(ni),lat(nj),time(ntime)
      real lono(nlon),lato(nlat)
      integer data_start(3),data_count(3)
      character*100 ref_date,time_calendar

      open(1, file='pca.cdl')
      write(1,*)'netcdf pca_loadings {'
      write(1,*)'dimensions:'
      write(1,*)'	factor = UNLIMITED ;'
      write(1,*)'	lon = ', nlon,' ;'
      write(1,*)'	lat = ', nlat,' ;'
      write(1,*)'       time = ', ntime,' ;'
      write(1,*)'       rclass = ', nclass,' ;'
      write(1,*)'       nwarea = ', nwarea,' ;'

      write(1,*)'variables:'
      write(1,*)'double factor(factor) ;'

      write(1,*)'double lon(lon) ;'
      write(1,*)'lon:units = "degrees_east" ;'
      write(1,*)'lon:long_name = "Lon" ;'

      write(1,*)'double lat(lat) ;'
      write(1,*)'lat:units = "degrees_north" ;'
      write(1,*)'lat:long_name = "Lat" ;'

      write(1,*)'double rclass(rclass) ;'
      write(1,*)'rclass:units =  "mm" ;'
      write(1,*)'rclass:long_name = "Rainfall intesity" ;'

      write(1,*)'double rfreq(rclass) ;'
      write(1,*)'rfreq:units =  "unit" ;'
      write(1,*)'rfreq:long_name = "Frequency" ;'

      write(1,*)'double rcumfreq(rclass) ;'
      write(1,*)'rcumfreq:units = "unit" ; '
      write(1,*)'rcumfreq:long_name = "Cumulative Frequency" ;'

      write(1,*)'double time(time) ;'
      write(1,*)'time:units = "', trim(ref_date), '" ;'
      write(1,*)'time:long_name = "Time" ;'
      write(1,*)'time:calendar = "' , trim(time_calendar), '" ;'

      write(1,*)'double nwarea(nwarea) ;'
      write(1,*)'nwarea:units =  "percent" ;'
      write(1,*)'nwarea:long_name = "percentage of area of extreme event" ;'

      write(1,*)'float Mask(lat, lon) ;'
      write(1,*)'Mask:long_name = "mask" ;'
      write(1,*)'Mask:units = "" ;'
      write(1,*)'Mask:missing_value = -1.0e+30;'
      write(1,*)'Mask:_FillValue = -1.0e+30;'

      write(1,*)'float Extreme(lat, lon) ;'
      write(1,*)'Extreme:long_name = "Extreme Treshold" ;'
      write(1,*)'Extreme:units = "" ;'
      write(1,*)'Extreme:missing_value = -1.0e+30;'
      write(1,*)'Extreme:_FillValue = -1.0e+30;'

      write(1,*)'float ExtremeDays(time, nwarea) ;'
      write(1,*)'ExtremeDays:long_name = "Extreme event over WCape" ;'
      write(1,*)'ExtremeDays:units = "" ;'
      write(1,*)'ExtremeDays:missing_value = -1.0e+30;'
      write(1,*)'ExtremeDays:_FillValue = -1.0e+30;'

      write(1,*)'float SOMS_Nodes(factor, lat, lon) ;'
      write(1,*)'SOMS_Nodes:long_name = "SOMs Nodes" ;'
      write(1,*)'SOMS_Nodes:units = "" ;'
      write(1,*)'SOMS_Nodes:missing_value = -1.0e+30;'
      write(1,*)'SOMS_Nodes:_FillValue = -1.0e+30;'

      write(1,*)'float Freq(factor) ;'
      write(1,*)'Freq:long_name = "Frequency" ;'
      write(1,*)'Freq:units = "" ;'
      write(1,*)'Freq:missing_value = -1.0e+30;'
      write(1,*)'Freq:_FillValue = -1.0e+30;'

      write(1,*)'float NodeTime(time) ;'
      write(1,*)'NodeTime:long_name = "The Node" ;'
      write(1,*)'NodeTime:units = "" ;'
      write(1,*)'NodeTime:missing_value = -1.0e+30;'
      write(1,*)'NodeTime:_FillValue = -1.0e+30;'

      write(1,*)'float pr(time, lat, lon) ;'
      write(1,*)'pr:long_name = "rainfall" ;'
      write(1,*)'pr:units = "mb" ;'
      write(1,*)'pr:missing_value = -1.0e+30;'
      write(1,*)'pr:_FillValue = -1.0e+30;'

  
      write(1,*)'}'

      close(1)      
      write(command,110) trim(OutFile)
      write(*,*) command
      call system(command)
     
      write(command,111) trim(OutFile)
      write(*,*) command
      call system(command)

  110 format('rm ' a50)
  111 format('ncgen -o ' a50, '   pca.cdl')
    

!---------------------------------------------------------------
! Open output file
!--------------------------------------------------------------
       iret = nf_open(OutFile,NF_WRITE,ncido)
       if(iret .ne. NF_NOERR) call handle_error(iret)

!---- Write lat
       jo = 0
       do j= lat_start,  lat_stop 
       jo = jo+1
       lato(jo) = lat(j)
       enddo
       iret = nf_inq_varid (ncido, 'lat', ivarid)
       if (iret.ne.NF_NOERR) call handle_error (iret)
       iret = nf_put_var_real (ncido,ivarid,lato)
       if (iret.ne.NF_NOERR) call handle_error(iret)

!---- Write lon
       io = 0
       do i= lon_start,lon_stop 
       io = io+1
       lono(io) = lon(i)
       enddo    
       iret = nf_inq_varid (ncido, 'lon', ivarid)
       if (iret.ne.NF_NOERR) call handle_error (iret)
       iret = nf_put_var_real (ncido,ivarid,lono)
       if (iret.ne.NF_NOERR) call handle_error(iret)


!---- Write lon
       iret = nf_inq_varid (ncido, 'time', ivarid)
       if (iret.ne.NF_NOERR) call handle_error (iret)
       iret = nf_put_var_real (ncido,ivarid,time)
       if (iret.ne.NF_NOERR) call handle_error(iret)



!---------------------------------------------------------------
! Close output file
!--------------------------------------------------------------

       iret = nf_close(ncido)
       if (ret.ne.NF_NOERR) then
       call handle_error (iret)
       endif

      end


!==================================================================
      Subroutine Prepare_CombinedSOMsOutputFile(OutFile, nlon, nlat, &
                   lon_start,  lon_stop, lat_start,  lat_stop, &
                   lon, lat, ni, nj, time, ntime, ndata,ref_date, time_calendar)

      include 'netcdf.inc'
      character*40 OutFile
      character*100 command
      real lon(ni),lat(nj)
      real lono(nlon),lato(nlat)
      integer data_start(3),data_count(3)
      character*100 ref_date, time_calendar


      open(1, file='pca.cdl')
      write(1,*)'netcdf pca_loadings {'
      write(1,*)'dimensions:'
      write(1,*)'	factor = UNLIMITED ;'
      write(1,*)'	lon = ', nlon,' ;'
      write(1,*)'	lat = ', nlat,' ;'
      write(1,*)'       time = ', ntime,' ;'
      write(1,*)'       dset = ', ndata,' ;'

      write(1,*)'variables:'
      write(1,*)'double factor(factor) ;'

      write(1,*)'double lon(lon) ;'
      write(1,*)'lon:units = "degrees_east" ;'
      write(1,*)'lon:long_name = "Lon" ;'

      write(1,*)'double lat(lat) ;'
      write(1,*)'lat:units = "degrees_north" ;'
      write(1,*)'lat:long_name = "Lat" ;'

      write(1,*)'double time(time) ;'
      write(1,*)'time:units = "', trim(ref_date), '" ;'
      write(1,*)'time:long_name = "Time" ;'
      write(1,*)'time:calendar = "' , trim(time_calendar), '" ;'


      write(1,*)'float SOMS_Nodes(factor, lat, lon) ;'
      write(1,*)'SOMS_Nodes:long_name = "SOMs Nodes" ;'
      write(1,*)'SOMS_Nodes:units = "" ;'
      write(1,*)'SOMS_Nodes:missing_value = -1.0e+30;'
      write(1,*)'SOMS_Nodes:_FillValue = -1.0e+30;'

      write(1,*)'float AllFreq(factor) ;'
      write(1,*)'AllFreq:long_name = "Frequency" ;'
      write(1,*)'AllFreq:units = "" ;'
      write(1,*)'AllFreq:missing_value = -1.0e+30;'
      write(1,*)'AllFreq:_FillValue = -1.0e+30;'

      write(1,*)'float SubFreq(factor,dset) ;'
      write(1,*)'SubFreq:long_name = "Frequency" ;'
      write(1,*)'SubFreq:units = "" ;'
      write(1,*)'SubFreq:missing_value = -1.0e+30;'
      write(1,*)'SubFreq:_FillValue = -1.0e+30;'

      write(1,*)'float NodeTime(time,dset) ;'
      write(1,*)'NodeTime:long_name = "The NodeTime" ;'
      write(1,*)'NodeTime:units = "" ;'
      write(1,*)'NodeTime:missing_value = -1.0e+30;'
      write(1,*)'NodeTime:_FillValue = -1.0e+30;'

      write(1,*)'}'

      close(1)      
      write(command,110) trim(OutFile)
      write(*,*) command
      call system(command)
     
      write(command,111) trim(OutFile)
      write(*,*) command
      call system(command)

  110 format('rm ' a50)
  111 format('ncgen -o ' a50, '   pca.cdl')
    

!---------------------------------------------------------------
! Open output file
!--------------------------------------------------------------
       iret = nf_open(OutFile,NF_WRITE,ncido)
       if(iret .ne. NF_NOERR) call handle_error(iret)

!---- Write lat
       jo = 0
       do j= lat_start,  lat_stop 
       jo = jo+1
       lato(jo) = lat(j)
       enddo
       iret = nf_inq_varid (ncido, 'lat', ivarid)
       if (iret.ne.NF_NOERR) call handle_error (iret)
       iret = nf_put_var_real (ncido,ivarid,lato)
       if (iret.ne.NF_NOERR) call handle_error(iret)

!---- Write lon
       io = 0
       do i= lon_start,lon_stop 
       io = io+1
       lono(io) = lon(i)
       enddo    
       iret = nf_inq_varid (ncido, 'lon', ivarid)
       if (iret.ne.NF_NOERR) call handle_error (iret)
       iret = nf_put_var_real (ncido,ivarid,lono)
       if (iret.ne.NF_NOERR) call handle_error(iret)

!---- Write lon
       iret = nf_inq_varid (ncido, 'time', ivarid)
       if (iret.ne.NF_NOERR) call handle_error (iret)
       iret = nf_put_var_real (ncido,ivarid,time)
       if (iret.ne.NF_NOERR) call handle_error(iret)



!---------------------------------------------------------------
! Close output file
!--------------------------------------------------------------

       iret = nf_close(ncido)
       if (ret.ne.NF_NOERR) then
       call handle_error (iret)
       endif

      end

