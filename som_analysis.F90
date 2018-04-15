        Program SOM
        
!       ifort -o som.exe -heap-arrays 1000000000 som_analysis.F90 -I/opt/netcdf/include/netcdf.inc -L/opt/netcdf/lib -lnetcdff
        include 'netcdf.inc'
!        parameter(ni=180,nj=88, ntime=10958,nfile=1,nvar=1) 
        parameter(ni=200,nj=140, ntime=7300,nfile=1,nvar=1) 
        character*13 vname(nvar)
        real lon(ni),lat(nj),time(ntime)
        character*40 OutFile1, OutFile
        character*40 lonName,latName,varName
        character*13 tname !vname(nvar),
        character*40 TimeUnitString
        character*100 command
        character*60 datafile
        character*5 TheLon, TheLat
        character*40 InFile(nfile)
        character*100 ref_date,time_calendar
!        real time(ntime)
        dimension ilev(nfile)
        real var(ni,nj)
        real rlonD_start, rlonD_stop,rlatD_start, rlatD_stop
        integer lonD_start,  lonD_stop, latD_start,  latD_stop
        integer nlonD, nlatD,nvar, idata, nxdim,nydim
        
        include "IO_specification.F90"
        !------Open Input file
        write(*,*) InFile(1)
        write(*,*) OutFile
        write(*,*) lonName
        write(*,*) latName
        write(*,*) tname
        write(*,*) varName
        iret = nf_open(InFile(1),NF_NOWRITE,ncid1)
        if(iret .ne. NF_NOERR) call handle_error(iret)
        call ReadLonLatTime(lonName,latName,lon,lat,ni,nj,ncid1)

        !------Close the file
        iret = nf_close(ncid1)
        if (ret.ne.NF_NOERR) then
        call handle_error (iret)
        endif
        write(*,*) (lon(i),i=1,10)
        write(*,*) (lat(i),i=1,10)

        !------------------------------------------------------------
        !     Set the study domain (i.e. Southern Africa)
        !------------------------------------------------------------
        call Set_Domain_Ids(lon, lat, ni, nj, &
        rlonD_start, rlonD_stop,rlatD_start, rlatD_stop, &
        lonD_start,  lonD_stop, latD_start,  latD_stop )
        nlonD = lonD_stop - lonD_start +1
        nlatD = latD_stop - latD_start +1
        
        write(*,*) "Gridsize of area of interest"
        write(*,*) "nlon=",nlonD
        write(*,*) "nlat=",nlatD
!        !------------------------------------------------
!        !    Prepare and open the output file (NC)
!        !-------------------------------------------------
        call Prepare_OutputFile(OutFile, nlonD, nlatD, &
        lonD_start,  lonD_stop, latD_start,  latD_stop, &
        lon, lat, ni, nj,time, ntime, ref_date, time_calendar)

        write(*,*)"Done creating ", trim(OutFile)


!        !----------------------------------------------
!        !     Read and write climate variable under
!        !     investigation 
!        !----------------------------------------------
!         lon1=lonD_start
!         lat1=latD_start
!ReadWriteAllVars(OutFile, InFile,ilev, vname, nlon, nlat, lon1, lat1, itime_start, ntime)
        ilev=0
        vname(1)=varName
!        do ivar = 1, nvar
!        write(*,*) vname(ivar)
!        Call ReadWriteAllVars(OutFile, InFile(ivar), ilev(ivar), vname(ivar), nlonD, nlatD, &
!        lonD_start, latD_start, itime_start, ntime)
!        write(*,*) "done for ", OutFile, ":", vname(ivar), nvar
!        enddo

!        !---------------------------------
!        !      Perform SOMs analysis
!        !---------------------------------
        idata=1
        itime_start=1
        call GenSOMs(OutFile, InFile(1),ilev(1), vname(1), nlonD, nlatD, lonD_start, &   
        latD_start, itime_start, ntime, idata,nxdim,nydim)
!GenSOMs(OutFile,InFile,ilev,vname, ni,nj,lon1,lat1,itime_start,ntime,idata)
        write(*,*) "Analysis successfully completed !!!"

        end program SOM

!======================================================================
!======================================================================
      Subroutine GenSOMs(OutFile,InFile,ilev,vname,&
      ni,nj,lon1,lat1,itime_start,ntime,idata,nxdim,nydim)
      include 'netcdf.inc'
      character*100 command
      character*40 OutFile
      character*40 InFile
      character*13 vname
      real var(ni,nj) !,ExtrDays(nwarea,ntime)
      character*15 TheDate(ntime)
      dimension pcvar(ni*nj), iDays(ntime)


        open(1, file= "x.asc")
        idt = 0
        nprow = 0

        open(8, file="iused_idx.dat")

!------Open Input file and read data
          iret = nf_open(InFile,NF_NOWRITE,ncid1)
          if(iret .ne. NF_NOERR) call handle_error(iret)
!        iwarea = 6
        do  itime = 1,ntime           
!         if (ExtrDays(iwarea,itime) .eq. 1) then
           idt = idt + 1
           itread = itime + (itime_start-1)
            call Read_RAWData(vname,var,ni, nj, ilev, lon1, lat1, itread, ncid1)
            if (trim(vname) .eq. 'psl') var = var/100.
            if (trim(vname) .eq. 'PRECIP') then
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
!               write(*,*) var(i,j)
             enddo
            enddo
!            write(*,*) (pcvar(i),i=10,80)
!            stop
!            write(*,*) "idt:",  idt, npcol
            if(idt .eq. 1) then        
              write(*,*) "npcol :", npcol  
              write(1,*) npcol   ! first line in SOM input file
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
!          endif 
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
!      nxdim = 3
!      nydim = 4
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
!====================================================================
     subroutine ReadWriteAllVars(OutFile, InFile,ilev, vname, nlon, nlat, lon1, lat1, itime_start, ntime)
      include 'netcdf.inc'
      real var(nlon,nlat)
      character*13 vname
      character*40 OutFile
      character*40 InFile
      integer data3_start(3),data3_count(3)


!------Open Input file and read data
       iret = nf_open(InFile,NF_NOWRITE,ncid1)
       if(iret .ne. NF_NOERR) call handle_error(iret)

!------Open output file
       iret = nf_open(OutFile,NF_WRITE,ncido)
       if(iret .ne. NF_NOERR) call handle_error(iret)

       do  itime = 1,ntime 
        itread = itime + (itime_start-1)
        call Read_RAWData(vname,var,nlon, nlat, ilev, lon1, lat1, itread, ncid1)
        if (trim(vname) .eq. 'psl') var = var/100.
        if (trim(vname) .eq. 'PR_FEWS')  var = var*86400.0

        data3_start(1) = 1
        data3_start(2) = 1
        data3_start(3) = itime

        data3_count(1) = nlon
        data3_count(2) = nlat
        data3_count(3) = 1 

        iret = nf_inq_varid (ncido, 'pr', ivarid)
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_put_vara_real (ncido,ivarid,data3_start,data3_count,var)
        if (iret.ne.NF_NOERR) call handle_error (iret)
       enddo

!===== close file
        iret = nf_close(ncid1)
        if (ret.ne.NF_NOERR) then
        call handle_error (iret)
        endif

        iret = nf_close(ncido)
        if (ret.ne.NF_NOERR) then
        call handle_error (iret)
        endif

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

!==================================================================================
      Subroutine ReadLonLatTime(lonName,latName,lon,lat,ni,nj,ncid1)
      include 'netcdf.inc'
      real lon(ni),lat(nj)!,time(ntime)
      real lon2D(ni,nj), lat2D(ni,nj)
      character*40 lonName,latName

!---- Read lat
        iret = nf_inq_varid (ncid1, trim(latName), ivarid)  
        if (iret.ne.NF_NOERR) call handle_error (iret)
        iret = nf_get_var_real (ncid1,ivarid,lat)
        if (iret.ne.NF_NOERR) call handle_error(iret)
!        lat = lat2D(1,:)

!---- Read lon
       iret = nf_inq_varid (ncid1, trim(lonName), ivarid)  
       if (iret.ne.NF_NOERR) call handle_error (iret)
       iret = nf_get_var_real (ncid1,ivarid,lon)
       if (iret.ne.NF_NOERR) call handle_error(iret)
!       lon = lon2D(:,1)       

!---- Read time
!       do it=1,ntime
!        itstart = it + (itime_start-1)
!        itcount = 1
!        iret = nf_inq_varid (ncid1, 'time', ivarid)  
!        if (iret.ne.NF_NOERR) call handle_error (iret)
!        iret = nf_get_vara_real (ncid1,ivarid,itstart,itcount,timein)
!        if (iret.ne.NF_NOERR) call handle_error(iret)
!        time(it) = timein
!        write(*,*) time(it)
!       enddo

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
                   lon, lat, ni, nj,time, ntime, ref_date, time_calendar)

      include 'netcdf.inc'
      character*40 OutFile
      character*100 command
      real lon(ni),lat(nj),time(ntime)
      real lono(nlon),lato(nlat)
      integer data_start(3),data_count(3)
      character*100 ref_date,time_calendar

      open(1, file='som.cdl')
      write(1,*)'netcdf som_loadings {'
      write(1,*)'dimensions:'
      write(1,*)'	factor = UNLIMITED ;'
      write(1,*)'	lon = ', nlon,' ;'
      write(1,*)'	lat = ', nlat,' ;'
      write(1,*)'       time = ', ntime,' ;'


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
  111 format('ncgen -o ' a50, '   som.cdl')
    

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

      end Subroutine Prepare_OutputFile

!===========================================================================
      subroutine handle_error(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
      print *, nf_strerror(iret)
      stop
      endif
      end

