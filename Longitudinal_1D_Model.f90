#include "AMYCAI.h"

MODULE Longitudinal_1D_Model

    USE TimeCatcher
    USE FloatPrecision
    USE ArraySize
    USE GRIDS
    USE hydro
    USE Parameters
    USE LinearSolver
    USE PRISM_RELA    
    
    IMPLICIT NONE

    contains

!= = = = = = = = = = = = = = = = = = = = = =
!
!   Main Routine to iterate the surface level
!   AND discharge
!
!= = = = = = = = = = = = = = = = = = = = = =  
    
    SUBROUTINE Saint_Venent_Model(T0, UU1, UU2, EE1, EE2, QQ1, ND, NDQ, &
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
    
        IMPLICIT NONE
        
        REAL(fp)                                          , INTENT(IN)    :: T0       
        REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: UU1
        REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: UU2
        REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: EE1 
        REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: EE2   
        REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: QQ1
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(IN)    :: BC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: HR
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET  , INTENT(IN)    :: FC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET  , INTENT(INOUT) :: HETA
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET  , INTENT(INOUT) :: Q_RESIDUAL2
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: DEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(INOUT) :: TDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(IN)    :: MDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(IN)    :: MAREA
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: MDEPE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: MTDEPE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: EtaU
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: MDEPU
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: MTDEPU
        INTEGER  , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: ND
        INTEGER  , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: NDQ
        
        REAL(fp) ,               DIMENSION(:) , POINTER                 :: TEMP
        LOGICAL                                                         :: GSTB
        INTEGER                                                         :: I
            
!- - - - - - - - - - - - - - - - - - - - - 
!
!   LPX: Solve the saint venent equation
!        AND control the flood AND dry
!
!- - - - - - - - - - - - - - - - - - - - - 
        LPX:DO
        
    ! Saint Venent equation
        
            CALL SaintVenent_Solver(T0,UU1,UU2,EE1,EE2,BC,HR,MTDEPU,MTDEPE,MAREA,FC,ND,NDQ)

    ! flood AND dry

            CALL VaryBnd(GSTB,EE1,EE2,MDEPE,MTDEPU,ND,NDQ)
        
    ! Keep iterating until there are no negative water depth

            IF(GSTB) THEN
                EXIT LPX
            END IF

        END DO LPX    
        
!- - - - - - - - - - - - - - - - - - - - -  
!   Update the mean water depth for next iteration
!- - - - - - - - - - - - - - - - - - - - -
        
        CALL Calculate_Total_Depth(UU2,EE2,DEP0,TDEP0,MDEPE,MDEPU,MTDEPE,MTDEPU,EtaU,BC,HR)

!- - - - - - - - - - - - - - - - - - - - -  
!   Calculate the cross-sectional discharge
!- - - - - - - - - - - - - - - - - - - - -
        
        !QQ1(:) = UU2(:)*MTDEPU(:)*BT
        
!- - - - - - - - - - - - - - - - - - - - -
! Calculate Maximum eta , Q_RESIDUAL
!- - - - - - - - - - - - - - - - - - - - -   
        
!        DO I = ML , MU
!            HETA(I)  = dmax1(HETA(I),EE2(I))
!#ifdef RESIDUALINC
!            Q_RESIDUAL2(I) = Q_RESIDUAL2(I) + QQ2(I)
!#endif
!        END DO

! - - - - - - - -
! Exchange Datas
! - - - - - - - -
        
        TEMP => EE1
        EE1  => EE2 
        EE2  => TEMP
        TEMP => UU1
        UU1  => UU2
        UU2  => TEMP  
        
! - - - - - - - -    
! release the temporal POINTER
! - - - - - - - -        
        NULLIFY(TEMP)
       
    END SUBROUTINE Saint_Venent_Model

!= = = = = = = = = = = = = = = = = = = = = =
! 
!          Saint Venent Equation
!
!= = = = = = = = = = = = = = = = = = = = = =  

    SUBROUTINE SaintVenent_Solver(T0,UU1,UU2,EE1,EE2,BC,HR,MTDEPU,MTDEPE,MAREA,FC,ND,NDQ)

        IMPLICIT NONE
    
        REAL(fp)                                        , INTENT(IN)    :: T0       
        REAL(fp)               , DIMENSION(:) , POINTER , INTENT(IN)    :: UU1
        REAL(fp)               , DIMENSION(:) , POINTER , INTENT(INOUT) :: UU2
        REAL(fp)               , DIMENSION(:) , POINTER , INTENT(IN)    :: EE1 
        REAL(fp)               , DIMENSION(:) , POINTER , INTENT(INOUT) :: EE2 
        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(IN)    :: BC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(IN)    :: HR
        REAL(fp) , ALLOCATABLE , DIMENSION(:) , TARGET  , INTENT(IN)    :: FC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(IN)    :: MTDEPU
        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(IN)    :: MTDEPE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(IN)    :: MAREA
        INTEGER  , ALLOCATABLE , DIMENSION(:)           , INTENT(IN)    :: ND        
        INTEGER  , ALLOCATABLE , DIMENSION(:)           , INTENT(IN)    :: NDQ   
        
        REAL(fp)                                                        :: UP 
        REAL(fp)                                                        :: ADVE , ADVW 
        REAL(fp)                                                        :: RESIST , DIFFX , CHEZY
        INTEGER                                                         :: I , J , K
        INTEGER                                                         :: JC , JJC , J0

!- - - - - - - - - - - - - - - - - - - - - 
!
!       CONTINUITY EQUATION LOOP: LPJ0
!
!- - - - - - - - - - - - - - - - - - - - - 
    
        JC = -1
        LPJ0:DO J = ML , MU+1 , 1
            
    ! ID in amt array
            
            JC = JC + 2
    
    ! boundary conditions
            
            IF(J==ML.OR.J==MU+1) THEN
                
    ! OPEN boundary: water level
                
                IF(ND(J)<0) THEN
                    AMT(1,JC) = 0.D0
                    AMT(2,JC) = 0.D0
                    AMT(3,JC) = 1.D0
                    AMT(4,JC) = 0.D0
                    AMT(5,JC) = 0.D0
                    AMT(6,JC) = BndSeaLevel(T0)
                    CYCLE LPJ0
                    
    ! closed boundary: no variation               
                    
                ELSEIF(ND(J)>0) THEN
                    AMT(1,JC) = 0.D0
                    AMT(2,JC) = 0.D0
                    AMT(3,JC) = 1.D0
                    AMT(4,JC) = 0.D0
                    AMT(5,JC) = 0.D0
                    AMT(6,JC) = EE1(J) 
                    CYCLE LPJ0   
                ENDIF       
    
    ! inner grids
    
            ELSE
                
    ! dry
                
                IF(ND(J)>0) THEN
                    AMT(1,JC) = 0.D0
                    AMT(2,JC) = 0.D0
                    AMT(3,JC) = 1.D0
                    AMT(4,JC) = 0.D0
                    AMT(5,JC) = 0.D0
                    AMT(6,JC) = EE1(J) 
                    CYCLE LPJ0
    ! wet
                    
                ELSE                   
                    
                    ! Consider the width change
                    
                    AMT(1,JC) =   0.D0
                    AMT(2,JC) = - DT*MTDEPU(J-1)    *BC(J-1)
                    AMT(3,JC) = + DY                * 0.5D0*(BC(J-1)+BC(J))               
                    AMT(4,JC) = + DT*MTDEPU(J)      *BC(J)
                    AMT(5,JC) =   0.D0
                    AMT(6,JC) = + DY*EE1(J)         * 0.5D0*(BC(J-1)+BC(J))
                    CYCLE LPJ0                   
                    
                ENDIF
            ENDIF  
        END DO LPJ0
        
!- - - - - - - - - - - - - - - - - - - - - 
!    
!        MOMENTUM EQUATION LOOP: LPJ1
!    
!- - - - - - - - - - - - - - - - - - - - - 
    
        JJC = 0
        LPJ1:DO J = ML , MU , 1

    ! ID in amt array
            
            JJC = JJC + 2
        
    ! boundary condition

            IF(J==ML.OR.J==MU) THEN
                
    ! closed bnd
                
                IF(NDQ(J)>0) THEN
                    AMT(1,JJC) = 0.D0
                    AMT(2,JJC) = 0.D0
                    AMT(3,JJC) = 1.D0
                    AMT(4,JJC) = 0.D0
                    AMT(5,JJC) = 0.D0
                    AMT(6,JJC) = 0.D0
                    CYCLE LPJ1
                    
    ! opened bnd
                    
                ELSEIF(NDQ(J)<=0) THEN
                    AMT(1,JJC) = 0.D0
                    AMT(2,JJC) = - GRA*DT
                    AMT(3,JJC) = + DY
                    AMT(4,JJC) = + GRA*DT
                    AMT(5,JJC) = 0.D0
                    AMT(6,JJC) = + DY*UU1(J)
                    CYCLE LPJ1     
                ENDIF     
                
            ENDIF          

    ! inner grids : dry case

            IF(NDQ(J)>0.OR.ND(J)>0.OR.ND(J+1)>0) THEN
                AMT(1,JJC) = 0.D0
                AMT(2,JJC) = 0.D0
                AMT(3,JJC) = 1.D0
                AMT(4,JJC) = 0.D0
                AMT(5,JJC) = 0.D0
                AMT(6,JJC) = 0.D0
                CYCLE LPJ1
            ENDIF
            
    ! inner grids: wet case (general case)
            
    ! Advection term

            ADVE = DMIN1(0.D0,UU1(J))
            ADVW = DMAX1(0.D0,UU1(J))
            
    ! Resistance term
            
            !RESIST = DABS(UU1(J))*GRA/65.D0**2/MTDEPU(J)
            CHEZY = 20.D0*HR(J)**(1.D0/6.D0)
            RESIST = DABS(UU1(J))*GRA/CHEZY**2/HR(J)
            !RESIST = FC(J)
            
    ! X-Diffusion term

            DIFFX = 50.D0

    ! Fill the matrix    
        
            AMT(1,JJC) = - ADVW*DT &
                       & - DIFFX*DT/DY
        
            AMT(2,JJC) = - GRA*DT

            AMT(3,JJC) = + DY &
                       & + (ADVW-ADVE)*DT &
                       & + RESIST*DT*DY &
                       & + 2.D0*DIFFX*DT/DY

            AMT(4,JJC) = + GRA*DT

            AMT(5,JJC) = + ADVE*DT &
                       & - DIFFX*DT/DY
        
            AMT(6,JJC) = + DY*UU1(J)
        
        END DO LPJ1

        JC = MAX0(JC,JJC)
    !    IF(I==59) THEN
    !OPEN(5,file = 'amt3.txt')
    !do J = 1 , jc
    !WRITE(5,'(20e20.10)') (amt(k,J),k=1,6)
    !end do
    !CLOSE(5)
    !WRITE(*,*) I
    !    END IF
        CALL SOLVE5(JC)
    
    !IF(I==59) THEN
    !OPEN(5,FILE = 'AMT4.TXT')
    !DO J = 1 , JC , 1
    !WRITE(5,'(20E20.10)') (AMT(K,J),K=1,6)
    !END DO
    !CLOSE(5)
    !PAUSE 'X'
    !    END IF

    ! farm
        
        J0 = 0
        DO J = 2 , JC , 2
            J0 = J0 + 1
            UU2(J0) = AMT(6,J)
        END DO
        J0 = 0
        DO J = 1 , JC , 2
            J0 = J0 + 1
            EE2(J0) = AMT(6,J)
        END DO
    
    END SUBROUTINE SaintVenent_Solver

!= = = = = = = = = = = = = = = = = = = = = =
! 
!          Flood AND Dry routine
!
!= = = = = = = = = = = = = = = = = = = = = =  

    SUBROUTINE VaryBnd(GSTB, EE1, EE2, MDEPE, MTDEPU, ND, NDQ)
    
        USE TimeCatcher
        USE FloatPrecision
        USE ArraySize
        USE Parameters

        IMPLICIT NONE

        LOGICAL                                         , INTENT(OUT)   :: GSTB      
        REAL(fp)               , DIMENSION(:) , POINTER , INTENT(IN)    :: EE1 
        REAL(fp)               , DIMENSION(:) , POINTER , INTENT(INOUT) :: EE2 
        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(INOUT) :: MDEPE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(INOUT) :: MTDEPU        
        INTEGER  , ALLOCATABLE , DIMENSION(:)           , INTENT(INOUT) :: ND        
        INTEGER  , ALLOCATABLE , DIMENSION(:)           , INTENT(INOUT) :: NDQ
        

        REAL(fp) :: H0 , S0
        INTEGER  :: I , CNODES , II , I1 , I2

        GSTB = .true.
        CNODES = 0

    ! Check IF there are negative water Depth (1D)
    
        DO I = ML+1 , MU , 1
            H0 = MDEPE(I) + EE2(I)
            IF(H0<=0.D0) THEN
                ND(I) = 1
                GSTB = .FALSE.
                CNODES = CNODES+1
            END IF
        END DO

        IF(.NOT.GSTB) RETURN

        DO I = ML+1 , MU , 1
            IF(ND(I)==1) THEN
                ND(I) = 0
            END IF
        END DO

    !- - - - - - - - - - - - - - - - - - - - 
    ! CHECK NDQ the status of velocity grids
    ! wet process: dry -> flood -> wet
    ! dry process: wet -> dry
    ! dry: 1
    ! flood: -5
    ! wet: 0
    !- - - - - - - - - - - - - - - - - - - - 
        
        DO I = ML , MU , 1
            
        ! wet grids 
        
            IF(NDQ(I)==0) THEN
                
        ! wet -> dry
                
                IF(MTDEPU(I)<=DDRY)  THEN
                    NDQ(I) = 1
                ENDIF
                
        ! flood grids
                
            ELSEIF(NDQ(I)==-5) THEN
                
        ! flood -> wet
                
                IF(MTDEPU(I)>=DWET)  THEN
                    NDQ(I) = 0
                ELSE             
                    
        ! flood -> dry
                    
                    S0 = (EE1(I)-EE1(I+1))*(EE2(I)-EE2(I+1))
                    IF(S0<0.D0) THEN
                        NDQ(I) = 1
                    ENDIF    
                    
                ENDIF
        
        ! dry grids
                
            ELSEIF(NDQ(I)==1) THEN
        
        ! dry -> wet
                
                IF(MTDEPU(I)>=DWET) THEN
                    NDQ(I) = 0
        
        ! dry -> flood
                    
                ELSE
                    IF(EE1(I)/=EE1(I+1)) THEN
                        S0 = (EE1(I)-EE1(I+1))*(EE2(I)-EE2(I+1))
                        IF(S0<0.D0) THEN
                            NDQ(I) = -5
                        ENDIF               
                    ELSE
                        IF(EE2(I)/=EE2(I+1)) THEN
                            NDQ(I) = -5
                        ENDIF         
                    ENDIF                    
                ENDIF
                
            ENDIF
            
        END DO

    END SUBROUTINE VaryBnd
   
!= = = = = = = = = = = = = = = = = = = = = =
! 
!   Calculate the Meandepth at eta grids and U grids
!
!= = = = = = = = = = = = = = = = = = = = = =  

    SUBROUTINE Calculate_Total_Depth(UU2,EE2,DEP0,TDEP0,MDEPE,MDEPU,MTDEPE,MTDEPU,EtaU,BC,HR)

        IMPLICIT NONE
    
        REAL(fp)               , DIMENSION(:)   , POINTER , INTENT(IN)    :: UU2 
        REAL(fp)               , DIMENSION(:)   , POINTER , INTENT(IN)    :: EE2     
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: DEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(INOUT) :: TDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(IN)    :: MDEPE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: MTDEPE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(IN)    :: MDEPU
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: MTDEPU
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: EtaU
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(IN) :: BC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: HR

        INTEGER             :: I , J
        INTEGER             :: status
        CHARACTER(len=80)   :: err_msg

! Total Dep at surface node    
    
        DO I = ML , MU+1
            MTDEPE(I) = MDEPE(I) + EE2(I)
        END DO

! Total Dep at velocity node

        DO I = ML , MU
            IF(UU2(I)>0.D0) THEN
                EtaU(I) = EE2(I)
            ELSEIF(UU2(I)<0.D0) THEN
                EtaU(I) = EE2(I+1)
            ELSE
                EtaU(I) = DMAX1(EE2(I),EE2(I+1))
            ENDIF
            MTDEPU(I) = MDEPU(I) + EtaU(I)    ! mean total dep
            TDEP0(I,:) = DEP0(I,:)   + EtaU(I)    ! local total dep
            HR(I) = MTDEPU(I)*BC(I)/(2.D0*MTDEPU(I)+BC(I))  ! Hydraulic Radius
        ENDDO
        
    END SUBROUTINE Calculate_Total_Depth
    
    SUBROUTINE Calculate_Steady_Depth(UU2, MDEP0, MDEPE, MDEPU)

        IMPLICIT NONE
    
        REAL(fp)               , DIMENSION(:) , POINTER , INTENT(INOUT) :: UU2
        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(INOUT) :: MDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(INOUT) :: MDEPE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(INOUT) :: MDEPU
        
        INTEGER                               :: I , J
        INTEGER                               :: status
        CHARACTER(len=80)                     :: err_msg

! Dep at surface node (MAX)   

        DO I = ML+1 , MU
            MDEPE(I) = DMAX1(MDEP0(I-1),MDEP0(I))
        END DO
        MDEPE(ML) = MDEPE(ML+1)
        MDEPE(MU+1) = MDEPE(MU)

! Dep at velocity node (UPWIND)
    
        DO I = ML , MU
            IF(UU2(I)>0.D0) THEN
                MDEPU(I) = MDEPE(I)
            ELSEIF(UU2(I)<0.D0) THEN
                MDEPU(I) = MDEPE(I+1)
            ELSE
                MDEPU(I) = DMIN1(MDEPE(I),MDEPE(I+1))
            ENDIF
        ENDDO
    
    END SUBROUTINE Calculate_Steady_Depth
    
!= = = = = = = = = = = = = = = = = = = = = =
! 
! Main Routine to solve highly simplified continuity equation
! to get the discharge and surface information
!
!= = = = = = = = = = = = = = = = = = = = = =  
    
    SUBROUTINE Water_Continuity_Model(EE1,DEE1,DEP0,T0,TDEP0,QQ1,BC,AC)
    
        IMPLICIT NONE
        
        REAL(fp)               , DIMENSION(:)   , POINTER , INTENT(INOUT) :: EE1           
        REAL(fp)               , DIMENSION(:)   , POINTER , INTENT(INOUT) :: DEE1
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(IN)    :: DEP0
        REAL(fp)                                          , INTENT(IN)    :: T0 
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(INOUT) :: TDEP0
        REAL(fp)               , DIMENSION(:)   , POINTER , INTENT(INOUT) :: QQ1 
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: BC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: AC

        REAL(fp)                          :: MSL , DMSL , PlanArea
        INTEGER                           :: I

#ifdef TIDE_OSCILLATION
    ! surface level
        DO I = ML , MU
            CALL  SeaLevel(T0,DEP0(I,NU-10),MSL,DMSL)
            EE1(I) = MSL
            DEE1(I) = DMSL
        END DO            
        !DMSL   = DBndSeaLevel(T0)          ! time derivative of the surface elevation       
        Call Calculate_Depth_Area_B(EE1,DEP0,TDEP0,BC,AC)   ! get the cross-sectional area and width of the water surface
!- - - - - - - - - - - -         
! compute the discharge
!- - - - - - - - - - - -          
        !PlanArea = 4.D0*BT**2
        PlanArea = 0.D0
        QQ1(MU) = PlanArea*DEE1(MU)
        DO I = MU-1 , ML , -1
            PlanArea = 0.5D0*(BC(I)+BC(I+1))*DY
            QQ1(I) = QQ1(I+1) + PlanArea*DEE1(I)
        END DO 

!tidal prism
#ifdef SSC_RELA_PA
        IF(T0>=0.5D0*TDCYC.AND.T0<=TDCYC) THEN
            TDPRISM = TDPRISM + QQ1(ML)*DT
        ENDIF
#endif

#elif defined(TIDE_MAX_DISCHARGE)

        EE1(:) = 0.D0 
        Call Calculate_Depth_Area_B(EE1,DEP0,TDEP0,BC,AC)  ! get the cross-sectional area and width of the water surface
        QQ1(MU) = Q_CONST
        DO I = MU-1 , ML , -1
            PlanArea = 0.5D0*(BC(I)+BC(I+1))*DY
            QQ1(I) = QQ1(I+1) + PlanArea*TDAMP*TDFRQ
        END DO     

#else

        EE1(:) = 0.D0               ! surface level always at mean water level        
        Call Calculate_Depth_Area_B(EE1,DEP0,TDEP0,BC,AC)  ! get the cross-sectional area and width of the water surface
        DO  I  = ML , MU
            QQ1(I) = DBLE(MU-I)/DBLE(MU-ML)*Q_CONST
        END DO  

#endif

    END SUBROUTINE Water_Continuity_Model

!- - - - - - - - - - - - - - -
!
! Boundary Condition Function
!
!- - - - - - - - - - - - - - - 

    SUBROUTINE SeaLevel(T0,D0,MSL,DMSL)

        IMPLICIT NONE
        REAL(fp)                                          , INTENT(IN)    :: D0
        REAL(fp)                                          , INTENT(IN)    :: T0     
        REAL(fp)                                          , INTENT(INOUT) :: MSL    
        REAL(fp)                                          , INTENT(INOUT) :: DMSL    
        REAL(fp)                                                          :: A1, B1, D1

#ifdef DYNAMIC_AMP
        D1 = D0
        A1 = 0.5D0*DMAX1(0.D0,(TDAMP+DMIN1(TDAMP,D1)))
        B1 = DMIN1(TDAMP,(A1 - DMIN1(D1,TDAMP)))
        MSL = B1 + A1*DCOS(TDFRQ*T0)
        DMSL = -A1*TDFRQ*DSIN(TDFRQ*T0)
#else
        MSL = TDAMP*DCOS(TDFRQ*T0)
        DMSL = -TDAMP*TDFRQ*DSIN(TDFRQ*T0)
#endif

    END SUBROUTINE SeaLevel
    
! - - - - - - - - - - - - - - 
    
    !function DBndSeaLevel(T0,D0)
    !
    !    IMPLICIT NONE
    !    REAL(fp)                                          , INTENT(IN)    :: D0
    !    REAL(fp)                                          , INTENT(IN)    :: T0         
    !    REAL(fp)                                                          :: A1, B1, D1, DBndSeaLevel
    !    REAL(fp) :: DBndSeaLevel, T0
    !
    !    DBndSeaLevel = -TDAMP*TDFRQ*DSIN(TDFRQ*T0)
    !
    !END function DBndSeaLevel    
    function BndSeaLevel(T0)
    
        IMPLICIT NONE
        REAL(fp)                                          , INTENT(IN)    :: T0         
        REAL(fp)                                                          :: BndSeaLevel
    
        BndSeaLevel = -TDAMP*DCOS(TDFRQ*T0)
    
    END function BndSeaLevel     
    
! - - - - - 
! calculate the area and width of the channel (continuity model)
! - - - - - 
    
    SUBROUTINE Calculate_Depth_Area_B(EE1,DEP0,TDEP0,BC,AC) 

        IMPLICIT NONE
        REAL(fp)               , DIMENSION(:)   , POINTER , INTENT(IN)    :: EE1           
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(IN)    :: DEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(INOUT) :: TDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: BC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: AC

        REAL(fp)                              :: CC , h1 , h2
        REAL(fp)                              :: aa1 , bb1
        INTEGER                               :: I , J
        INTEGER                               :: status
        CHARACTER(len=80)                     :: err_msg
            
! separate method
        
! calculate the area and width
        
        DO I = ML , MU

! Total Water Depth
            
            TDEP0(I,:) = DEP0(I,:) + EE1(I)
                    
! Width and Area
            
            aa1 = 0.d0        
            bb1 = 0.d0    
            do J = nl , nu-1
                h1 = TDEP0(I,J)
                h2 = TDEP0(I,J+1)
#ifndef TIDE_OSCILLATION
                h1 = h1 - 0.005D0
                h2 = h2 - 0.005D0
#endif
                if(h1>0.d0.and.h2>0.d0) then
                    aa1 = aa1 + 0.5d0*dx*(h1+h2)
                    bb1 = bb1 + dx
                elseif(h1>0.d0.and.h2<=0.d0) then
                    CC = h1 / (DEP0(I,J) - DEP0(I,J+1))
                    aa1 = aa1 + 0.5d0*dx*h1*CC
                    bb1 = bb1 + dx*CC
                elseif(h1<=0.d0.and.h2>0.d0) then
                    CC = h2 / (DEP0(I,J+1) - DEP0(I,J))
                    aa1 = aa1 + 0.5d0*dx*h2*CC
                    bb1 = bb1 + dx*CC
                end if                        
            end do           
            AC(I) = aa1
            BC(I) = bb1
        end do
          
    END SUBROUTINE Calculate_Depth_Area_B
    
!---------
    
    SUBROUTINE Calculate_TDAREA(DEP0,TDAREA) 

        IMPLICIT NONE
        REAL(fp)                                                          :: TDAREA        
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(IN)    :: DEP0

        REAL(fp)                              :: CC
        REAL(fp)                              :: aa1 , bb1
        INTEGER                               :: I , J
        INTEGER                               :: status
        CHARACTER(len=80)                     :: err_msg
            
! separate method
        
! calculate the area and width
        
        DO I = ML , ML

! Width and Area
            
            aa1 = 0.d0            
            do J = nl , nu-1
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
            TDAREA = aa1
        end do
          
    END SUBROUTINE Calculate_TDAREA        
    
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

END MODULE Longitudinal_1D_Model


    
    