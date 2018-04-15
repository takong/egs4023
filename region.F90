!==== Set and select region of interest
! This file is used to select the region of interest for the computations

   if (iregion .eq. 1) then
!===================================
!   For Southern Africa (Drakensberg)
!===================================
       ! Southern Africa
       rlonD_start =  0
       rlonD_stop  =  50
       rlatD_start = -45
       rlatD_stop  = -10

       !Drakensberg
       rlonA_start =  26   ! start lon
       rlonA_stop  =  32    ! end lon
       rlatA_start =  -28   ! start lat
       rlatA_stop  =  -32   ! end lat

       OutFile1  =  "dmr_output.nc"
       iwcmask  =  1
       ivar_used    = 1
       trs_per_extr = 60

    else
       ! Southern Africa
       rlonD_start =  0
       rlonD_stop  =  50
       rlatD_start = -45
       rlatD_stop  = -10

       !Drakensberg
       rlonA_start =  26   ! start lon
       rlonA_stop  =  32    ! end lon
       rlatA_start =  -28   ! start lat
       rlatA_stop  =  -32   ! end lat

       OutFile1  =  "other_output.nc"
       iwcmask  =  1
       ivar_used    = 1
       trs_per_extr = 60
        
    endif



