!
!  Reservoir scheme of LEAF-Hydro-Flood-Dam (LHFD) model
!
!  Shin, S., Y. Pokhrel, and G. Miguez-Macho (2019) 
!     “High Resolution Modeling of Reservoir Release and Storage Dynamics at the Continental Scale”, 
!     Water Resources Research, Vol. 55, Issue 1, pp.787-810. doi:10.1029/2018WR023025
!
!  Originally authored by Gonzalo Miguez-Macho (https://scholar.google.com/citations?user=wSKBQLsAAAAJ).
!  Modified by Sanghoon Shin (https://sites.google.com/view/sshin) & Yadu Pokhrel (https://water.egr.msu.edu/)
!
!
MODULE module_leafhydro
    use module_parallel ! subroutines for MPI programing

    implicit none

    integer, parameter :: nDAMs=2250
    integer, save :: RESs, RESe
	
! nDAMs			total number of dams [-]
! RESs			index of dams to be simulated (start) [-]
! RESe			index of dams to be simulated (end) [-]

CONTAINS

! NOTE:

! At every routing timestep (dtlr=60sec)
! 	call LHFD_routing
!
! At every days,
! 	call update_damSTOR
!
! At every months,
! 	call update_k_rls
! 	call update_prov_r
! 	call update_plan_r_versatile
	
subroutine LHFD_routing(imax,js,je,frwtd,dtll,dtlr,veg,input_fd,input_bfd,qnew,qs,qrf,qspring &
                ,sriver,slope,depth,width,length,maxdepth,area,topo_rivers,floodheight,flood,delsfcwat,qmean &
                ,damID,damYEAR,damCap,damSRF,damH,damPRP,damFLAG,damMINSTOR,damRv,resID,damSTOR,table_STOR,downdiff &
                ,plan_r)
				
! imax			x-coordinate [-]
! js			y-coordinate (start) [-]
! je			y-coordinate (end) [-]
! frwtd			timestep of ground water modules [s]
! dtll			timestep of land surface modules [s]
! dtlr			timestep of routing module [s]
! veg			not used
! input_fd		flow direction towards downstream [-]
! input_bfd		flow direction towards upstream [-]
! qnew			river discharge [m3/s]
! qs			direct runoff [m3/m2]
! qrf			groundwater inflow to river [m3/m2]  
! qspring		spring from groundwater to surface [m3/m2]
! sriver		river water depth [m]
! slope			riverbed slope [-]
! depth			river water depth [m]
! width			river width [m]
! length		river length [m]
! maxdepth		river depth to start brim over [m]
! area			area of grid cell [m2]
! topo_rivers	river water elevation with bankful discharge [m]
! floodheight	flooded water depth [m]
! flood			flooded water flow [m3/s]
! delsfcwat		snow water exchange [m]
! qmean			mean river discharge (for diagnose purpose)
! damID			dam ID [-]
! damYEAR		dam comission year [-]
! damCap		dam storage capacity [m3]
! damSRF		dam surface area [km2]
! damH			dam height [m]
! damPRP		dam purpose [-]
! damFLAG		not used
! damMINSTOR	the minimum reservoir storage [%]
! damRv			dam release coefficient [-]
! resID			reservoir ID [-]
! damSTOR		dam storage [m3]
! table_STOR	dam storage table [m3]
! downdiff		river bed slope [-]
! plan_r		target release (rm) [m3/s]

    implicit none

    integer, parameter :: ntsplit=2
    integer :: i,j,imax,js,je,n,k,nveg,i1,j1,i2,j2,iter
    integer ,dimension(imax,js:je), intent(in) :: input_fd,input_bfd
    integer ,dimension(imax,js:je) :: fd,bfd
    real, dimension(imax,js:je) :: veg,q,qnew,qs,qrf,qspring,qext,velold &
        ,slope,depth,width,length,maxdepth,area,topo_rivers,floodheight,flood &
        ,waterelev,sriver,delsfcwat,qmean,plan_r
    real :: dtll,frwtd,dtlr,deltat,   aa,wi,speed,depthmax,dflood,deltadepth,areaflood,areariver,dsnew &
        ,dswstor,watslope,hcommon,factor,aaa,bbb,ccc,dx,v1,v2,v
    integer :: reqsu,reqsd,reqru,reqrd
    real(8) ::         swstorold(imax,js:je),riverchannel
    real(8) :: swstor
    real :: qin,surf_area_m3,suplus_stor

    integer, dimension(imax,js:je), intent(in) :: damID, &      
                                                damYEAR, &      
                                                damPRP, &       
                                                damFLAG, &      
                                                resID            
    real, dimension(imax,js:je), intent(in):: damCap, &         
                                            damSRF, &           
                                            damH, &             
                                            damMINSTOR, &       
                                            damRv               
    real, dimension(imax,js:je), intent(in out):: damSTOR       
    real, dimension(imax,js:je), intent(in) :: downdiff         
    real tempstorage
    real, dimension(nDAMs,RESs:RESe), intent(in out) :: table_STOR

    real, dimension(imax,js:je) ::  areariver_array,&   
                                    areaflood_array     
    real, dimension(imax,js:je) :: alpha, &                       !
                                 spill                          !
    real :: relCoeff, &             
            minCap, &               
            CurrentStorage          
    
    
    fd = input_fd
    bfd = input_bfd
    


    
    !input of water in the cell            
    qext = ( (qrf+qspring)/frwtd + (qs+delsfcwat)/dtll ) * area + flood



    do n=1,ntsplit 


        do iter=1,2

IF((pid.gt.0.and.pid.lt.numtasks-1).or.numtasks.eq.1)then 
            if(numtasks.gt.1)call sendborders(imax,js,je,qnew,reqsu,reqsd,reqru,reqrd)

            !make sure that the borders are received before calculating anything
            if(pid.eq.1)then
                call  MPI_wait(reqru,status,ierr)
            elseif(pid.eq.numtasks-2)then
                call  MPI_wait(reqrd,status,ierr)
            elseif(pid.gt.1.and.pid.lt.numtasks-2)then
                call  MPI_wait(reqru,status,ierr)
                call  MPI_wait(reqrd,status,ierr)
            endif

            q=qnew



            if(iter.eq.1)then
                deltat=0.5*dtlr/float(ntsplit) 
                where(sriver.gt.0.)
                    velold=qnew/(width*sriver) 
                elsewhere
                    velold=0.
                endwhere
            else
                deltat=dtlr/float(ntsplit)
            endif


            do j=js+1,je-1
                do i=3,imax-2

                    IF(fd(i,j).ne.0 )then
                    !calculate total inflow into cell i j
                    !fd (flow direction) tells where the river in cell i j is flowing to

                        qin=0.

                        if(fd(i+1,j+1).eq.  8)qin=qin+q(i+1,j+1) !NE
                        if(fd(i  ,j+1).eq.  4)qin=qin+q(i  ,j+1) !N
                        if(fd(i-1,j+1).eq.  2)qin=qin+q(i-1,j+1) !NW
                        if(fd(i-1,j  ).eq.  1)qin=qin+q(i-1,j  ) !W
                        if(fd(i-1,j-1).eq.128)qin=qin+q(i-1,j-1) !SW
                        if(fd(i  ,j-1).eq. 64)qin=qin+q(i  ,j-1) !S
                        if(fd(i+1,j-1).eq. 32)qin=qin+q(i+1,j-1) !SE
                        if(fd(i+1,j  ).eq. 16)qin=qin+q(i+1,j  ) !E


                        dsnew = qin-q(i,j)  !SHIN: qin is sum of inflows from upstreams, q(i,j) is outflow to a downstream.

                        !old surface water storage
                        areariver=width(i,j)*length(i,j)
                        areaflood = max( area(i,j)-areariver, 0. )

                        if(iter.eq.1)then
                            swstorold(i,j) = dble(sriver(i,j))*dble(areariver) + dble(floodheight(i,j))*dble(areaflood)
                        endif

                        !new surface water storage
                        swstor = swstorold(i,j) + ( dble(dsnew)+dble(qext(i,j)) ) *dble(deltat)

                        !now redistribute between river channel and floodplain
                        riverchannel=dble(maxdepth(i,j))*dble(areariver)


                        if(swstor.gt.riverchannel)then
                            floodheight(i,j)=(swstor-riverchannel)/max(area(i,j),areariver)
                            sriver(i,j)=floodheight(i,j)+maxdepth(i,j)
                        else
                            floodheight(i,j)=0.
                            sriver(i,j) = swstor/areariver
                        endif
                    ENDIF
                enddo
            enddo


            !before changing qnew make sure that the borders have been received
            if(pid.eq.1)then
                call  MPI_wait(reqsu,status,ierr)
            elseif(pid.eq.numtasks-2)then
                call  MPI_wait(reqsd,status,ierr)
            elseif(pid.gt.1.and.pid.lt.numtasks-2)then
                call  MPI_wait(reqsu,status,ierr)
                call  MPI_wait(reqsd,status,ierr)
            endif


            if(numtasks.gt.1)call sendborders(imax,js,je,sriver,reqsu,reqsd,reqru,reqrd)

            !make sure that the borders are received before calculating anything
            if(pid.eq.1)then
                call  MPI_wait(reqru,status,ierr)
            elseif(pid.eq.numtasks-2)then
                call  MPI_wait(reqrd,status,ierr)
            elseif(pid.gt.1.and.pid.lt.numtasks-2)then
                call  MPI_wait(reqru,status,ierr)
                call  MPI_wait(reqrd,status,ierr)
            endif

            depth = max(sriver,0.)

            areariver_array=width*length                     !areariver=width(i,j)*length(i,j)
            areaflood_array=area-areariver_array             !areaflood = max( area(i,j)-areariver, 0. )
            where (areaflood_array < 0.) areaflood_array=0.

            ! if the bfd cell is dam, set bfd=0
            do j=js+1,je-1
                do i=3,imax-2
                    call flowdir(imax,js,je,bfd,i,j,i2,j2)
                    if (i2==0 .and. j2==0) cycle
                    if (damID(i2,j2)/=0) bfd(i,j)=0
                enddo
            enddo  
            
            do j=js+1,je-1
                do i=3,imax-2
					! For dams
                    if(damID(i,j)/=0)then
                    
                        if (Res_option == 0) then  ! Doll 2003
                              qnew(i,j) = 1./(24.*60.*60.) * 0.01 * damSTOR(i,j) * ( damSTOR(i,j)/damCap(i,j) )**1.5
                        elseif (Res_option >= 1 ) then  ! "1:Hanasaki" , "2:Biemans ", "3:Rold_M01", "4:Rnew_M01", "5:Rcal_M01"
                              qnew(i,j)=plan_r(i,j)
                        endif
                        
                        ! For all Res_option....
                        if (sriver(i,j)/damH(i,j) < 0.1) qnew(i,j) = 0.
                        ! (a2) the minimum release :: damSTOR should not below damMINSTOR
                        if (damSTOR(i,j) < damMINSTOR(i,j)) qnew(i,j) = 0.
                        ! (b) the spilling
                        if (damSTOR(i,j) > damCAP(i,j)) then !The surplus --> released for 1 day
                              qnew(i,j) = qnew(i,j) + 1./(24.*60.*60.) * (1./1.) * ( damSTOR(i,j) - damCAP(i,j) )
                        elseif (sriver(i,j) - damH(i,j) > 0.001) then
                              surf_area_m3 = max(damSRF(i,j),0.1) * 10.**6.  
                              suplus_stor = surf_area_m3 * (sriver(i,j)-damH(i,j))
                              qnew(i,j) = qnew(i,j) + 1./(24.*60.*60.) * (1./1.) * suplus_stor * damSTOR(i,j)/damCAP(i,j)
                        endif
                    else
                    !  For non-dams
                        IF(fd(i,j).gt.0 )then
                            !calculate residence time based on velocity for present flow
                            if(sriver(i,j).gt.1.e-9)then
                                call flowdir(imax,js,je,fd,i,j,i1,j1)
                                dx=0.5*(length(i,j)+length(i1,j1))
                                call flowdir(imax,js,je,bfd,i,j,i2,j2)

                                ! if the next cell is reservoir cell, adjust watslope
                                if ( resID(i,j) == 0 .and. resID(i1,j1) /= 0) then
                                    watslope=(                              - depth(i,j) ) /dx
                                else
                                    watslope=( downdiff(i,j) + depth(i1,j1) - depth(i,j) ) /dx
                                endif

                                if(j1.eq.1.or.j1.eq.n3big)then
                                    aa=depth(i,j)*width(i,j)/(2.*depth(i,j)+width(i,j))
                                    speed = ( aa**(2./3.) )*sqrt(slope(i,j))/0.03
                                    speed=min(max(speed,0.001),length(i,j)/deltat)
                                    qnew(i,j)=width(i,j)*depth(i,j)*speed
                                else
                                    if(fd(i1,j1).ne.0)then

                                        if(watslope.lt.0.)then
                                             aa=depth(i,j)
                                        else                   
                                             aa=max(depth(i1,j1)- downdiff(i,j),0.)
                                        endif

                                        if(aa.gt.0.)then

                                            if(velold(i,j).gt.0.)then
                                                if (i2.eq.0.or.j1.eq.0) then 
                                                    aaa=1.+deltat*( (velold(i,j))/length(i,j) + g*0.03**2. * velold(i,j)/aa**(4./3.) )
                                                    if ( j1 == 0) print *, "check SHIN.. it is possilble?", i,j,i1,j1,i2,j2, fd(i,j), bfd(i,j) 
                                                else 
                                                    aaa=1.+deltat*( (velold(i,j)-velold(i2,j2))/length(i,j) + g*0.03**2. * velold(i,j)/aa**(4./3.) )
                                                endif
                                                ccc=-deltat*g*watslope+velold(i,j)
                                            else
                                                aaa=1.+deltat*( (velold(i1,j1)-velold(i,j))/length(i1,j1) + g*0.03**2. * velold(i,j)/aa**(4./3.) )
                                                ccc=-deltat*g*watslope+velold(i,j)
                                            endif
                                        else
                                            aaa=1.
                                            ccc=0.
                                        endif

                                    else
                                        aa=depth(i,j)

                                        if(velold(i,j).gt.0.)then
                                            if (i2.eq.0.or.j1.eq.0) then
                                                aaa=1.+deltat*( (velold(i,j))/length(i,j) + g*0.03**2. * velold(i,j)/aa**(4./3.) )
                                            else 
                                                aaa=1.+deltat*( (velold(i,j)-velold(i2,j2))/length(i,j) + g*0.03**2. * velold(i,j)/aa**(4./3.) )
                                            endif
                                            ccc=-deltat*g*watslope+velold(i,j)
                                        else
                                            aaa=1.+deltat*(  g*0.03**2. * velold(i,j)/aa**(4./3.) )
                                            ccc=-deltat*g*watslope+velold(i,j)
                                        endif

                                    endif

                                    qnew(i,j)= width(i,j)*depth(i,j)*ccc/aaa

                                    if(qnew(i,j).gt.0.)then
                                        qnew(i,j)=min(qnew(i,j),width(i,j)*depth(i,j)*10.)
                                    elseif(fd(i1,j1).ne.0)then
                                        if(depth(i1,j1).gt.0.)then
                                            qnew(i,j)=-min(-qnew(i,j),width(i,j)*depth(i,j)*1.)
                                        else
                                            qnew(i,j)=0.
                                        endif
                                    endif
                                endif
                            else
                                qnew(i,j)=0.
                            endif
                            
                        else
                            qnew(i,j) = 0.
                        ENDIF
                    endif
                enddo
            enddo

            !accumulate river flow to do daily mean later !SHIN: qmeam is not used here though
            if(iter.eq.1)qmean=qmean+ ( qnew*dtlr/float(ntsplit) )

            !before changing sriver make sure that the borders have been received
            if(pid.eq.1)then
                call  MPI_wait(reqsu,status,ierr)
            elseif(pid.eq.numtasks-2)then
                call  MPI_wait(reqsd,status,ierr)
            elseif(pid.gt.1.and.pid.lt.numtasks-2)then
                call  MPI_wait(reqsu,status,ierr)
                call  MPI_wait(reqsd,status,ierr)
            endif


ENDIF 

        enddo
    enddo
end subroutine LHFD_routing


!******************************************************************************************
subroutine FLOWDIR(imax,js,je,fd,ii,jj,i,j)
    implicit none
    integer :: imax,js,je,i,j,ii,jj
    integer, dimension(imax,js:je) :: fd

    select case(fd(ii,jj))
    case(2,4,8)
        j=jj-1
    case(1,16)
        j=jj
    case(32,64,128)
        j=jj+1
    case default
        j=0
    end select

    select case(fd(ii,jj))
    case(128,1,2)
        i=ii+1
    case(4,64)
        i=ii
    case(8,16,32)
        i=ii-1
    case default
        i=0
    end select

end subroutine flowdir


subroutine update_damSTOR(imax,js,je,resID,width,length,area,sriver,floodheight,damID,table_STOR,damSTOR)
    integer :: i,j,imax,js,je
    
    integer, dimension(imax,js:je), intent(in) :: resID,damID  ! reservoir extent map .. corresponding damID of dam map
    real, dimension(imax,js:je), intent(in) :: width, length,area,sriver,floodheight
    real, dimension(nDAMs,RESs:RESe), intent(in out) :: table_STOR
    real, dimension(imax,js:je), intent(out) :: damSTOR

    real :: tempstorage
    real, dimension(imax,js:je) ::  areariver_array,&   
                                    areaflood_array     
                                      
    areariver_array=width*length                     
    areaflood_array=area-areariver_array             
    where (areaflood_array < 0.) areaflood_array=0.
    
    table_STOR = 0.
    damSTOR=0.       !reset WHOLE 'damSTOR' for every calls
    do j=js+1,je-1
       do i=3,imax-2
          if (resID(i,j)/= 0) then
             tempstorage = sriver(i,j)*areariver_array(i,j) + floodheight(i,j)*areaflood_array(i,j)
             table_STOR(resID(i,j),pid) = table_STOR(resID(i,j),pid) + tempstorage
          endif
       enddo
    enddo
    
    call Table_receive_and_sendback(table_STOR)
    
    do j=js+1,je-1
      do i=3,imax-2
        if(damID(i,j) /= 0) then
            damSTOR(i,j) = table_STOR(resID(i,j),pid)
        endif
      enddo
    enddo

end subroutine update_damSTOR

subroutine update_k_rls(n2,js,je,damID,damCAP, damSTOR,dam_OPM, monthforc,k_rls)
    integer             :: i,j,k
    real                :: alpha
    integer, intent(in) :: n2,js,je,monthforc
    real,    dimension(n2,js:je), intent(in) :: damCAP, damSTOR
    integer, dimension(n2,js:je), intent(in) :: dam_OPM,damID
    real,    dimension(n2,js:je),intent(inout) :: k_rls
    
!    alpha = 0.85
    alpha = 0.65
    
    do j=js,je
        do i=1,n2
            if (damID(i,j) /= 0) then
                if (dam_OPM(i,j) == monthforc) then
                    k_rls(i,j) = damSTOR(i,j) / (alpha*damCAP(i,j))
                endif
            endif
        enddo
    enddo
end subroutine update_k_rls

subroutine update_prov_r(n2,js,je,damID,DPI,dam_avgQ,SW_annual,SW_monthly,prov_r)
    integer             :: i,j,k
    real                :: alpha, M_eqsel, M_release
    integer, intent(in) :: n2,js,je
    integer, dimension(n2,js:je), intent(in) :: damID
    real,    dimension(n2,js:je), intent(in) :: DPI, dam_avgQ, SW_annual,SW_monthly
    real,    dimension(n2,js:je),intent(out) :: prov_r
    
    alpha = 0.85
    prov_r= 0.
    if (Res_option == 0) then      !0: Doll,
        M_eqsel = 0.5; M_release = 0.5
    elseif (Res_option == 1) then  !1:Hanasaki [M=0.5 DPI>0.5 R=max(1,(c/0.5)^2)],
        M_eqsel = 0.5; M_release = 0.5
    elseif (Res_option == 2) then  !2:Biemans [M=0.1 DPI>0.5 R=max(1,(c/0.5)^2)],
        M_eqsel = 0.5; M_release = 0.1
    elseif (Res_option == 3) then  !3:Rold_M01 [M=0.1 DPI>0.9 R=max(1,(c/0.5)^2)],
        M_eqsel = 0.9; M_release = 0.1
    elseif (Res_option == 4) then  !4:Rnew_M01 [M=0.1 DPI>0.9 R=max(1,(c/(1/0.85))^2)],
        M_eqsel = 0.9; M_release = 0.1
    elseif (Res_option == 5) then  !5:Rcal_M01 [M=0.1 DPI>0.9 R=max(1,(c/(1/0.85))^2) & calibrated with M=0.1],
        M_eqsel = 0.9; M_release = 0.1
    endif

    do j=js,je
        do i=1,n2
            if (damID(i,j) /= 0) then
                if (SW_annual(i,j) == 0.) then
                        prov_r(i,j) = abs(dam_avgQ(i,j))
                else
                    if (DPI(i,j) < M_eqsel) then
                        prov_r(i,j) = abs(dam_avgQ(i,j)) +SW_monthly(i,j) - SW_annual(i,j)
                    else
                        prov_r(i,j) = abs(dam_avgQ(i,j)) * (M_release + (1-M_release)*SW_monthly(i,j)/SW_annual(i,j) )
                    endif
                endif
            endif
        enddo
    enddo
        
end subroutine update_prov_r

subroutine update_plan_r_versatile(n2,js,je,damID,k_rls,dam_monQ,c_coeff,prov_r,monthforc,plan_r,damRv)
    integer             :: i,j,k
    real                :: c_threshold,alpha,const
    integer, intent(in) :: n2,js,je
    integer, intent(in) :: monthforc
    integer, dimension(n2,js:je), intent(in) :: damID
    real,    dimension(n2,js:je), intent(in) :: k_rls, c_coeff,prov_r
    real,    dimension(12,n2,js:je), intent(in) :: dam_monQ
    real,    dimension(n2,js:je),intent(in) :: damRv
    real,    dimension(n2,js:je),intent(out) :: plan_r
    
    
    alpha = 0.85
    plan_r= 0.
    if (Res_option == 0) then      !0:storage-based,
        !SHIN20171102 :: Now, it is OK
        !stop "STOPPED :: storage-based option (Res_option=0) is not available"
        c_threshold = 0.5      ; const = 2.    !=1/c_threshold   
    elseif (Res_option == 1) then  !1:Hanasaki [M=0.5 DPI>0.5 R=max(1,(c/0.5)^2)],
        c_threshold = 0.5      ; const = 2.    !=1/c_threshold
    elseif (Res_option == 2) then  !2:Biemans [M=0.1 DPI>0.5 R=max(1,(c/0.5)^2)],
        c_threshold = 0.5      ; const = 2.    !=1/c_threshold
    elseif (Res_option == 3) then  !3:Rold_M01 [M=0.1 DPI>0.9 R=max(1,(c/0.5)^2)],
        c_threshold = 1/alpha  ; const= alpha  !=1/c_threshold
    elseif (Res_option == 4) then  !4:Rnew_M01 [M=0.1 DPI>0.9 R=max(1,(c/(1/0.85))^2)],
        c_threshold = 1/alpha  ; const= alpha  !=1/c_threshold
    elseif (Res_option == 5) then  !5:Rcal_M01 [M=0.1 DPI>0.9 R=max(1,(c/(1/0.85))^2) & calibrated with M=0.1],
        c_threshold = 1/alpha  ; const= alpha  !=1/c_threshold
    endif
    
    do j=js,je
        do i=1,n2
            if (damID(i,j) /= 0) then
                if (Res_option <= 4) then !Non-Calibrated
                    if ( abs(c_coeff(i,j)) > c_threshold) then  !if c_coeff is negative, plan_r has huge magnitude.
                        plan_r(i,j) = k_rls(i,j) * prov_r(i,j)
                    else
                        plan_r(i,j) = abs(c_coeff(i,j))**2. * const**2. * k_rls(i,j) * prov_r(i,j) &
                                        + (1.- abs(c_coeff(i,j))**2. * const**2.)* abs(dam_monQ(monthforc,i,j))
                    endif
                else  ! Calibrated... here, use R-value calibrated to M=0.1
                    plan_r(i,j) = damRv(i,j) * k_rls(i,j) * prov_r(i,j) &
                                    + (1.- damRv(i,j))* abs(dam_monQ(monthforc,i,j))
                endif
            endif
        enddo
    enddo
        
end subroutine update_plan_r_versatile

