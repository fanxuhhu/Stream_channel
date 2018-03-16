! = = = = = = = = = = = = = = = = = = = = = = = 
!
!  This SUBROUTINE calculates the PA relationship
!
! = = = = = = = = = = = = = = = = = = = = = = = 

#include "AMYCAI.h"
    
SUBROUTINE QA_RELATION

! parameter AND array modules
  
    USE FloatPrecision
    USE TimeCatcher
    USE ArraySize
    USE Grids
    USE Hydro
    USE Topography
    USE Parameters
    
! Routine modules
    
    USE Longitudinal_1D_Model

    IMPLICIT NONE

! Define the Pointers to the surface AND 
! discharge array for 1D Longitudinal model
    
    REAL(fp) , DIMENSION(:)   , POINTER :: QQ1 , QQ2 
    REAL(fp) , DIMENSION(:)   , POINTER :: EE1 , EE2 , DEE1

    REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: TD_Prism 
    REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: Steady_Area
    
    REAL(fp)                          :: MAXH, aa1, CC
    REAL                              :: TT1, TT2, TTC
    INTEGER                           :: status , I , J , FID , FID0
    CHARACTER(len=100)                 :: FILENAME
    CHARACTER(len=80)                 :: DIRECTION
    CHARACTER(len=80)                 :: err_msg

    NULLIFY(QQ1,QQ2,EE1,EE2,DEE1)

! CPU Timer : begin
    
    CALL CPU_TIME(TT1)

! - - - - - - - - - - - 
!
! Initialize parameters
!
! - - - - - - - - - - -
    BT = DBLE(NU-NL)*dx
    LU = DBLE(MU-ML)*dy   
    T0 = 0.5D0*TDCYC   ! start from low tide
    IT0   = 0
    IHOUR = 0
    IDAY  = 0
    IYEAR = 0
    IDT   = INT(DT)
    CD    = 0
    CD2   = 0
    FID0 = 0
    
    Q_RESIDUAL1(:) = 0.D0
    Q_RESIDUAL2(:) = 0.D0
    
! ALLOCATE
    
    ALLOCATE(TD_Prism(ML:MU,1:File_Number),stat=status,errmsg=err_msg)
    ALLOCATE(Steady_Area(ML:MU,1:File_Number),stat=status,errmsg=err_msg)
    
! TARGET READY

    ETA1(:) = 0.D0
    ETA2(:) = 0.D0
    Q1(:) = 0.D0
    Q2(:) = 0.D0
    HETA(:) = 0.D0
    Q_RESIDUAL1(:) = 0.D0
    Q_RESIDUAL2(:) = 0.D0
    TDEP0 = 0.D0
    TD_Prism(:,:) = 0.D0
    Steady_Area(:,:) = 0.D0
    
!POINTER READY

    QQ1=>Q1
    QQ2=>Q2
    EE1=>ETA1
    EE2=>ETA2    
    DEE1=>DETA1
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = =       
!
!       Calculate the PA relationship
!
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    
! initialize the tidal amp
    
    TDAMP = TDAMP0

    DIRECTION = 'c:\tidal_channel_q2d_model_xf\Q2D_Model_2_longitudinal\XUFAN_QUASI_2D_model\'
    WRITE(*,*) File_Number
    pause
    LPF: DO FID = Start_File_ID , End_File_ID , File_ID_FRQ
        
        FID0 = FID0 + 1
        
        WRITE(*,*) 'FILE_ID = ' , FID , FID0
        WRITE(FILENAME,'(A76,A4,I6.6,A4)') TRIM(DIRECTION),'Dep_',FID,'.RES'
        OPEN(UNIT=10,FILE=FILENAME,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=status,ACTION='READ')
        READ(10) DEP0
        CLOSE(0)
        
!- - - - - - - - - - - - - - - - - - - - - - - - - - -         
!   Calculate the Area
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        DO I = ML , MU
            aa1 = 0.d0           
            do J = NL , NU-1
                if(DEP0(I,J)>0.d0.and.DEP0(I,J+1)>0.d0) then
                    aa1 = aa1 + 0.5d0*dx*(DEP0(I,J)+DEP0(I,J+1))
                elseif(DEP0(I,J)>0.d0.and.DEP0(I,J+1)<=0.d0) then
                    CC = DEP0(I,J) / (DEP0(I,J) - DEP0(I,J+1))
                    aa1 = aa1 + 0.5d0*dx*DEP0(I,J)*CC
                elseif(DEP0(I,J)<=0.d0.and.DEP0(I,J+1)>0.d0) then
                    CC = DEP0(I,J+1) / (DEP0(I,J+1) - DEP0(I,J))
                    aa1 = aa1 + 0.5d0*dx*DEP0(I,J+1)*CC
                end if                        
            end do           
            Steady_Area(I,FID0) = aa1*2.D0
#ifndef TIDE_OSCILLATION
            !Steady_Area(I,FID0) = Steady_Area(I,FID0) - 0.05D0*BT
#endif
        END DO       

        LPT:DO

!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!   Time variation
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
        
            T0 = T0 + DT
        
!- - - - - - - - - - - - - - - - - - - - - - - - - - -         
!   Longitudinal 1D model, calculate discharge Q 
!   , surface elevation eta AND cross-section area
!- - - - - - - - - - - - - - - - - - - - - - - - - - -     
#ifdef TIDE_OSCILLATION
            CALL Water_Continuity_Model(EE1,DEE1,DEP0,T0,TDEP0,QQ1,BC,AC)
#else
            DO  I  = ML , MU
                QQ1(I) = DBLE(MU-I)/DBLE(MU-ML)*Q_CONST
            END DO     
#endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - -         
!   Calculate the tidal prism
!- - - - - - - - - - - - - - - - - - - - - - - - - - -         

            TD_Prism(:,FID0) = TD_Prism(:,FID0) + QQ1(:)*dt*2.D0
            
!- - - - - - - - - - - - - - - - - - - - - - - - - - -   
! reaching the high tide , exit
!- - - - - - - - - - - - - - - - - - - - - - - - - - -         
        
            IF(T0>=TDCYC) THEN
                T0 = 0.5D0*TDCYC
                EXIT LPT
            END IF        

        END DO LPT
    END DO LPF
    
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

! OUTPUT
    
    OPEN(UNIT=10,FILE=TRIM(DIRECTION)//'Prism2.RES',FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=status,ACTION='WRITE')
    WRITE(10) TD_Prism
    CLOSE(10)

    OPEN(UNIT=10,FILE=TRIM(DIRECTION)//'Area2.RES',FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=status,ACTION='WRITE')
    WRITE(10) Steady_Area
    CLOSE(10)
    
    DEALLOCATE(TD_Prism,stat=status,errmsg=err_msg)
    DEALLOCATE(Steady_Area,stat=status,errmsg=err_msg)
    
! CPU Timer : END
    
    CALL CPU_TIME(TT2)
    WRITE(*,101) (TT2-TT1)/60.
    PAUSE

101 FORMAT("The simulation costs ",f5.1," minuts")
102 FORMAT("DAY = ",I5)
    
END SUBROUTINE QA_RELATION
