! = = = = = = = = = = = = = = = = = = = = = = = 
!
!  This SUBROUTINE is the kernal of the model
!  , which solves the 2D flow field, sediment
!  transport AND the bed evolution...
!
! = = = = = = = = = = = = = = = = = = = = = = = 

#include "AMYCAI.h"
    
SUBROUTINE Morphodynamic

! parameter AND array modules
  
    USE FloatPrecision
    USE TimeCatcher
    USE ArraySize
    USE Grids
    USE Hydro
    USE Topography
    USE Parameters
    USE ARTIFICIAL
    USE PRISM_RELA
    USE OUTPUT
    
! Routine modules
    USE Geology_Information 
    USE Longitudinal_1D_Model
    USE Lateral_1D_Model
    USE Sediment_Transport

    IMPLICIT NONE

! Define the Pointers to the surface AND 
! discharge array for 1D Longitudinal model


    REAL(fp) , DIMENSION(:)   , POINTER :: UU1 , UU2     
    REAL(fp) , DIMENSION(:)   , POINTER :: QQ1 , QQ2 
    REAL(fp) , DIMENSION(:)   , POINTER :: EE1 , EE2
    REAL(fp) , DIMENSION(:)   , POINTER :: DEE1
    REAL(fp) , DIMENSION(:)   , POINTER :: SSC1D1 , SSC1D2 
    REAL(fp) , DIMENSION(:,:) , POINTER :: SSC1 , SSC2 
        
    REAL(fp)                          :: MAXH
    REAL                              :: TT1, TT2, TTC
    INTEGER                           :: status , I , J
    CHARACTER(len=80)                 :: err_msg

    NULLIFY(UU1,UU2,QQ1,QQ2,EE1,EE2,DEE1)

! CPU Timer : begin
    
    CALL CPU_TIME(TT1)
    
! - - - - - - - - - - - 
!
! Initialization
!
! - - - - - - - - - - -

! Time Catcher
    
    T0 = 0.D0
    IT0   = 0
    IHOUR = 0
    IDAY  = 0
    IYEAR = 0
    IDT   = INT(DT)
    CD    = 0
    CD2   = 0

    Q_RESIDUAL1(:) = 0.D0
    Q_RESIDUAL2(:) = 0.D0

! TARGET READY
    
    DETA1(:) = 0.D0
    ETA1(:) = 0.D0
    ETA2(:) = 0.D0
    U1(:) = 0.D0
    U2(:) = 0.D0
    Q1(:) = 0.D0
    Q2(:) = 0.D0
    SC2D1(:,:) = 0.D0
    SC2D2(:,:) = 0.D0
    SC1D1(:) = 0.D0
    SC1D2(:) = 0.D0
    EVOSTEP(:,:) = 0.D0
    Evolution_TDCYC(:,:) = 0.D0
    HETA(:) = 0.D0
    HMTDEP0(:) = 0.D0
    MaxTau(:,:) = 0.D0
    Q_RESIDUAL1(:) = 0.D0
    Q_RESIDUAL2(:) = 0.D0
    TDEP0 = 0.D0

!POINTER READY

    UU1=>U1
    UU2=>U2
    QQ1=>Q1
    QQ2=>Q2
    EE1=>ETA1
    EE2=>ETA2
    DEE1=>DETA1
    SSC1=>SC2D1
    SSC2=>SC2D2
    SSC1D1=>SC1D1
    SSC1D2=>SC1D2

! - - - - - - - - - - - 
!
! Initialize morphodynamic information
!
! - - - - - - - - - - -

    
#ifdef ARTIFICIAL_EVO
    CALL ARTIFICIAL_BC_MDEP
    CALL Calculate_Steady_Depth(UU2, MDEP0, MDEPE, MDEPU)
    CALL Calculate_Total_Depth(UU2,EE2,DEP0,TDEP0,MDEPE,MDEPU,MTDEPE,MTDEPU,EtaU,BC,HR)
    CALL Grid_ND
    CALL Calculate_Total_Depth(UU2,EE2,DEP0,TDEP0,MDEPE,MDEPU,MTDEPE,MTDEPU,EtaU,BC,HR)
#elif defined(SAINT_VENENT)
    CALL Bathymetry
    CALL BedParameters( DEP0, MDEP0, DDEP, Upsilon, MAREA)
    CALL Sedimentation_Parameters(DEP0,TDEP0,HTDEP0,MDEP0,HMTDEP0,HETA,HND,HNDD,HDIFF)
    CALL Calculate_Steady_Depth(UU2, MDEP0, MDEPE, MDEPU)
    CALL Calculate_Total_Depth(UU2,EE2,DEP0,TDEP0,MDEPE,MDEPU,MTDEPE,MTDEPU,EtaU,BC,HR)
    CALL Grid_ND
    CALL Calculate_Total_Depth(UU2,EE2,DEP0,TDEP0,MDEPE,MDEPU,MTDEPE,MTDEPU,EtaU,BC,HR)
#elif defined(CONTINUITY)
    CALL Bathymetry
    CALL BedParameters( DEP0, MDEP0, DDEP, Upsilon, MAREA)
! PRISM
#ifdef SSC_RELA_PA
    TDPRISM0 = 2.D0*TDAMP0*BT*LU
    WRITE(*,*) TDAMP0,BT,LU
    pause
    PRISMRATIO = 1.D0
#endif
#endif 
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = =       
!
!       Flow model : 1D longitudinal AND 1D lateral
!
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    
! initialize the tidal amp
    
#ifdef INCREASE_AMP
    TDAMP = 0.002D0
#else
    TDAMP = TDAMP0
#endif
    
    LPT:DO

! Time variation
        
        T0 = T0 + DT

!- - - - - - - - - - - - - - - - - - - - - - - - - - -         
! Longitudinal 1D model, calculate discharge Q 
! , surface elevation eta AND cross-section area
!- - - - - - - - - - - - - - - - - - - - - - - - - - -     

#ifdef SAINT_VENENT  
        CALL Saint_Venent_Model(T0, UU1, UU2, EE1, EE2, QQ1, ND, NDQ, &
            &                   BC,                                     &            
            &                   HR,                                     &            
            &                   FC,                                     &
            &                   DEP0,                                   &   
            &                   TDEP0,                             &        
            &                   MDEP0,                               &
            &                   MAREA,                                       &
            &                   MDEPE,                             &
            &                   MTDEPE,                        &
            &                   EtaU,                                   &
            &                   MDEPU,                               &
            &                   MTDEPU,                          &
            &                   HETA,                                &
            &                   Q_RESIDUAL2                             &
            &                                                           )
#elif defined(CONTINUITY)
        CALL Water_Continuity_Model(EE1,DEE1,DEP0,T0,TDEP0,QQ1,BC,AC)
        CALL Sedimentation_Parameters(DEP0,TDEP0,HTDEP0,MDEP0,HMTDEP0,HETA,HND,HNDD,HDIFF)
#else
        WRITE(*,*) "Please Choose a flow model in file AMYCAI.h"                
#endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - -         
! Lateral 1D model, calculate transversal distribution 
! of velocity TDUD AND shear stress tau, AND get the
! corrected friction coefficient
!- - - - - - - - - - - - - - - - - - - - - - - - - - - 
        CALL Cross_Sectional_Tau(TDUD,TDTAU,QQ1,TDEP0,Upsilon,XDistance,WetSegments,FC,MaxTau,MeanTau,Sf,RHOGSD)    
        
#ifdef WAVE_INC
        CALL WAVE_INDUCED_TAU(WVTAU)
#endif
        
!- - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!           sediment transport model
!        
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifndef MAXTAUEVO
        CALL MUD_Transport(TDUD,TDTAU,MAXTAU,CRTAU,STDEP,HND,HNDD,HDIFF,HMTDEP0,Q_RESIDUAL1,SSC1,SSC2,SSC1D1,SSC1D2,ERO,TTERO,ERO1D,LABL,LAAVA,EVOSTEP,DEP0,UPSILON,BC,Sf,T0,TDEP0)
#endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!           change the seabed parameters
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - 
#ifndef MAXTAUEVO
#ifdef LONGTERM
        !IF(CD2 >= 15) THEN ! Active the evolution
#ifdef ARTIFICIAL_EVO
            CALL ARTIFICIAL_BC_MDEP_EVO(UU1)
#else
            CALL BedParameters( DEP0, MDEP0, DDEP, Upsilon, MAREA)
#endif
        !END IF
#endif  
#endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - -   
! After each tidal cycle, output
!- - - - - - - - - - - - - - - - - - - - - - - - - - -   
        IF(T0>=TDCYC) THEN
! Resetting time
            T0 = T0-TDCYC
            
! evolution
#ifdef MAXTAUEVO
#ifdef LONGTERM    
            CALL MUD_Transport(TDUD,TDTAU,MAXTAU,CRTAU,STDEP,HND,HNDD,HDIFF,HMTDEP0,Q_RESIDUAL1,SSC1D1,SSC1D2,ERO,TTERO,ERO1D,LABL,LAAVA,EVOSTEP,DEP0,UPSILON,BC,Sf,T0,TDEP0)
            CALL BedParameters( DEP0, MDEP0, DDEP, Upsilon, MAREA)
#endif
#endif
! Gradually increase the amplitude (linearly)     
#ifdef INCREASE_AMP
            if(TDAMP<TDAMP0) THEN
                TDAMP = TDAMP + 0.002D0
            endif
            if(TDAMP>TDAMP0) THEN
                TDAMP = TDAMP0
            endif
#endif
! change the SSC based on TIDAL PRISM
#ifdef SSC_RELA_PA
            Call Calculate_TDAREA(DEP0,TDAREA) 
            PRISMRATIO = (TDPRISM/TDPRISM0)!**(0.9D0)
            !PRISMRATIO = 1.576D-4*TDPRISM**0.95D0/TDAREA
            !IF (PRISMRATIO>1.D0) PRISMRATIO = 1.D0/PRISMRATIO
            TDPRISM = 0.D0
#endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - -     
! Update residual discharge
!- - - - - - - - - - - - - - - - - - - - - - - - - - -   
!#ifdef RESIDUALINC
!            Q_RESIDUAL1(:) = Q_RESIDUAL2(:)*DT/TDCYC
!            Q_RESIDUAL2(:) = 0.D0
!#endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - -            
! Output the topography data each tidal cycle
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
            CD2 = CD2 + 1
            WRITE(*,*) CD2 , dep0(1,nu),maxtau(1,nu)
            write(20,*) dep0(1,nu),maxtau(1,nu)
            !WRITE(*,*) DEP0(1,NL),DEP0(1,NU)
            !WRITE(50,'(4f20.10)') QQ1(1),QQ1(2),BC(1),BC(2)
            
#ifdef AVALANCHE2
            CALL AVALANCHE(DEP0)
            CALL BedParameters( DEP0, MDEP0, DDEP, Upsilon, MAREA)
#endif
            
#ifdef LONGTERM

#ifdef ARTIFICIAL_EVO
            BC(:) = BC(:) + BCEVO(:)
            WHERE(BC<0.1D0)
                BC = 0.1D0
            ENDWHERE
            !CALL Output_Binary_1D('DEP_','.RES',CD2,MDEP0)  
            CALL Output_Binary_1D('BCC_','.RES',CD2,BC)  
            BCEVO(:) = 0.D0
#else
            MeanTAU(:,:) = MeanTAU(:,:) / tdcyc * dt
            TTERO(:,:) = TTERO(:,:) / Morfactor
            CALL Output_Binary_2D('Dep_','.RES',CD2,DEP0)    
            !CALL Output_Binary_2D('LAB_','.RES',CD2,LABL)  
            CALL Output_Binary_2D('AAU_','.RES',CD2,MeanTau)
            !CALL Output_Binary_2D('MAU_','.RES',CD2,MaxTau)
            !CALL Output_Binary_2D_P('SSC_','.RES',CD2,SSC1)
            !CALL Output_Binary_2D('ERO_','.RES',CD2,TTERO)
            !CALL Output_Binary_2D('RGD_','.RES',CD2,RHOGSD)
            !CRTAUL(:,:) = CRTAU(1,:,:)
            !CALL OUTPUT_BInary_2D('CRT_','.RES',CD2,CRTAUL)  
            !WRITE(199,*) DEP0(1,NL),DEP0(1,NU)
            MaxTau(:,:)  = 0.D0
            MeanTau(:,:) = 0.D0
            TTERO(:,:) = 0.D0
#endif
            !IF(CD2>5000) EXIT LPT
#endif
        END IF        
        !WRITE(201,*) SSC1(4,NL)
!- - - - - - - - - - - - - - - - - - - - - - - - - - -      
! Output the instantaneous velocity , shear stress , SSC et.al 
!- - - - - - - - - - - - - - - - - - - - - - - - - - -        
             CD = CD + 1
#ifdef SHORTTERM        
            !CD = CD + 1
            !WRITE(*,*) CD            
! Output binary file for instantaneous variables
        IF(MOD(CD,1)==0) THEN
            CD3 = CD
            write(16,*) tdtau(1,nu-1)
            !CALL Output_Binary_1D('GAM_','.RES',CD,FC)
            !CALL Output_Binary_1D('ETA_','.RES',CD3,EE1)
            !CALL Output_Binary_1D('AAA_','.RES',CD,A1)
            !CALL Output_Binary_1D_INTEGER('NDD_','.RES',CD,NDQ)
            !CALL Output_Binary_1D('QQQ_','.RES',CD3,QQ1)
            !CALL Output_Binary_1D('UUU_','.RES',CD,SSC1D1)
            !CALL Output_Binary_1D('MTV_','.RES',CD,TDTAU(:,(NL+NU)/2))
            !CALL Output_Binary_2D_P('SSC_','.RES',CD3,SSC1)
            !CALL Output_Binary_1D('SS3_','.RES',CD,SSC1D1)
            !CALL Output_Binary_2D('DDD_','.RES',CD,TDTAU)
            !CALL Output_Binary_2D('ERO_','.RES',CD,ERO)
            !CALL Output_Binary_2D('EVO_','.RES',CD,EVOSTEP)
            !CALL Output_Binary_2D('HHH_','.RES',CD,DEP0) 

            
        ENDIF
! 2000 steps
            !WRITE(89,*) TDTAU(1,101)
            IF(CD>2000) EXIT LPT
#endif
  

!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!        
!   Simulation Timer: control the period AND output
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - -        

! update time

        IT0 = IT0 + IDT       
        
! T > 3600s
        
        IF(IT0>=SECONDS) THEN
            
            IT0 = IT0-SECONDS
            IHOUR = IHOUR + 1
            !WRITE(*,*) IHOUR
       
            
! T > 24h
            IF(IHOUR>=HOURS) THEN
                
                IHOUR = IHOUR - HOURS
                IDAY  = IDAY + Int(MorFactor)

!- - - - - - - - - - 
!  Output
!- - - - - - - - - - 
                
                !CD2 = CD2 + 1
                !WRITE(*,102) CD2
                !CALL Output_Binary_2D('Dep_','.RES',CD2,DEP0)
                !IF(IDAY>3) exit LPT

! T > 360days
                
                IF(IDAY>=DAYS) THEN
                    
                    IDAY = IDAY - DAYS
                    IYEAR = IYEAR + 1

                    
! T reach the upper limit years, terminate the simulation            
                    
                    IF(IYEAR>=YEARS+1) THEN
                        
                        exit LPT
                        
                    END IF
                    
                END IF
                
            END IF
            
        END IF
        
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
        
    END DO LPT

!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    
! CPU Timer : END
    
    CALL CPU_TIME(TT2)
    WRITE(*,101) (TT2-TT1)/60.
    !PAUSE

101 FORMAT("The simulation costs ",f15.5," Minutes")
102 FORMAT("DAY = ",I5)
    
    
END SUBROUTINE Morphodynamic
