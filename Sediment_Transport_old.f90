! This MODULE contains all the routines that 
! are relative to sediment transport. Only the diffusion
! proccess is considered for the 1st verson.    

    
#include "AMYCAI.h"
    
MODULE Sediment_Transport

    USE TimeCatcher
    USE FloatPrecision
    USE ArraySize
    USE GRIDS
    USE hydro
    USE Topography
    USE Parameters
    USE Suspended_Load
    USE LinearSolver
    
    IMPLICIT NONE

    REAL(fp) , DIMENSION(:) , POINTER :: VARS, VARP, VARN
    REAL(fp) , DIMENSION(:) , POINTER :: VSCS, VSCP, VSCN
    INTEGER  , DIMENSION(:) , POINTER :: GIDS, GIDP, GIDN
    REAL(fp) , DIMENSION(:) , POINTER :: CEP , HP  , FRUIT
    REAL(fp) , DIMENSION(:) , POINTER :: REQP 
    REAL(fp) , POINTER                :: LX , LY  , LT
    INTEGER  , POINTER                :: ILOW , IUP , JLOW , JUP
    REAL(fp)                          :: LXX , LYY

    REAL(fp) , DIMENSION(:) , POINTER :: BP
    
    CONTAINS


! = = = = = = = = = = = = = = = =    
!
!      1D diffusion model
!
! = = = = = = = = = = = = = = = =

    SUBROUTINE SSC_Diffusion_1D_Model(SSC1_1D,SSC2_1D)
    
        IMPLICIT NONE
        REAL(fp) , DIMENSION(:) , POINTER :: SSC1_1D , SSC2_1D , TEMP_1D
    
! Calculate local erosion        
    
        CALL Local_Erosion
        CALL Local_Total_Erosion
        
! 1D diffusion equation

! Y-Direcsion
 
        LX   => DY
        LT   => DT
        LXX = LX**2
        
        ILOW => ML
        IUP  => MU
        JLOW => NL
        JUP  => NU
            
! results of last step
            
        VARP  => SSC1_1D(:)
            
! wet AND dry information
            
        GIDP  => NDD_HWL_1D(:)
            
! diffusion coefficient -Y
            
        VSCP  => SSC_DIFFU_HWL_Y_1D(:)
            
! erosion term: local source term
            
        CEP   => ERO_1D(:)
            
! local water Depth
            
        HP    => Max_Area(:)

! width

        BP    => B_HWL(:)

! Residual Discharge

        REQP  => Q_RESIDUAL1(:)

! result
            
        FRUIT => SSC2_1D(:)     

! 1D solver
            
        CALL Tridiagonal_1D('Y',1)

! Exchange pointers    
    
        TEMP_1D=>SSC2_1D
        SSC2_1D=>SSC1_1D
        SSC1_1D=>TEMP_1D
        NULLIFY(VARP,FRUIT,GIDP,VSCP,TEMP_1D,CEP,BP,REQP) 
        NULLIFY(LX,LT,ILOW,IUP)
        
    END SUBROUTINE SSC_Diffusion_1D_Model
        
! = = = = = = = = = = = = = = = = 
!
!  Diffusion proccess for SSC
!
! = = = = = = = = = = = = = = = =     
   
    SUBROUTINE SSC_Diffusion_2D_Model(SSC1,SSC2)

        IMPLICIT NONE

        REAL(fp) , DIMENSION(:,:) , POINTER :: SSC1 , SSC2 , TEMP
        INTEGER :: i , j , i1 , i2 , j1 , j2


! Calculate local erosion

        CALL Local_Erosion

!- - - - - - - - - - - - - - 
!    ADI method is used
!- - - - - - - - - - - - - - 
    
! Y-Direcsion
 
        LX   => DY
        LY   => DX
        LT   => DT
        LXX = LX**2
        LYY = LY**2
        
        ILOW => ML
        IUP  => MU
        JLOW => NL
        JUP  => NU
        
        DO j = JLOW , JUP
            
            j1 = j+1
            j2 = j-1
            IF(j==JLOW) j2=j
            IF(j==JUP) j1=j
            
! results of last step
            
            VARS  => SSC1(:,j2)
            VARP  => SSC1(:,j)
            VARN  => SSC1(:,j1)
            
! wet AND dry information
            
            GIDS  => NDD_HWL(:,j2)
            GIDP  => NDD_HWL(:,j)
            GIDN  => NDD_HWL(:,j1)
            
! diffusion coefficient -Y
            
            VSCP  => SSC_DIFFU_HWL_Y(:,j)
            
! diffusion coefficient -X
            
            VSCS  => SSC_DIFFU_HWL_X(:,j2)
            VSCN  => SSC_DIFFU_HWL_X(:,j)
            
! erosion term: local source term
            
            CEP   => ERO(:,j)
            
! local water Depth
            
            HP    => Dep_HWL(:,j)

! result
            
            FRUIT => SSC2(:,j)     

! ADI solver
            
            CALL Tridiagonal_ADI('Y',j)
            
        END DO

! Exchange pointers    
    
        TEMP=>SSC2
        SSC2=>SSC1
        SSC1=>TEMP
        NULLIFY(VARS,VARP,VARN,FRUIT,GIDS,GIDP,GIDN,VSCP,VSCS,VSCN,TEMP,CEP) 
        NULLIFY(LX,LY,LT,ILOW,IUP,JLOW,JUP)
 
 ! - - - - - - - - - - - - - - - - -

 ! X-Direcsion

        LX   => DX
        LY   => DY
        LT   => DT
        LXX = LX**2
        LYY = LY**2
        
        ILOW => NL
        IUP  => NU
        JLOW => ML
        JUP  => MU
        
        DO j = JLOW , JUP
            
            j1 = j+1
            j2 = j-1
            IF(j==JLOW) j2=j
            IF(j==JUP ) j1=j
            
! results of last step
            
            VARS  => SSC1(j2,:)
            VARP  => SSC1(j,:)
            VARN  => SSC1(j1,:)
            
! wet AND dry information
            
            GIDS  => NDD_HWL(j2,:)
            GIDP  => NDD_HWL(j,:)
            GIDN  => NDD_HWL(j1,:)
            
! diffusion coefficient -Y
            
            VSCP  => SSC_DIFFU_HWL_X(j,:)
            
! diffusion coefficient -X
            
            VSCS  => SSC_DIFFU_HWL_Y(j2,:)
            VSCN  => SSC_DIFFU_HWL_Y(j,:)
            
! erosion term: local source term
            
            CEP   => ERO(j,:)
            
! local water Depth
            
            HP    => Dep_HWL(j,:)

! result
            
            FRUIT => SSC2(j,:)     

! ADI solver
            
            CALL Tridiagonal_ADI('X',j)
            
        END DO

! Exchange pointers    
    
        TEMP=>SSC2
        SSC2=>SSC1
        SSC1=>TEMP
        NULLIFY(VARS,VARP,VARN,FRUIT,GIDS,GIDP,GIDN,VSCP,VSCS,VSCN,TEMP,CEP) 
        NULLIFY(LX,LY,LT,ILOW,IUP,JLOW,JUP)

    END SUBROUTINE SSC_Diffusion_2D_Model
    
! = = = = = = = = = = = = = = = =    
!
!     SSC constant model
!
! = = = = = = = = = = = = = = = =  
    
    SUBROUTINE SSC_Constant_Model(SSC1)
        IMPLICIT NONE

        REAL(fp) , DIMENSION(:,:) , POINTER :: SSC1
        
        CALL Local_Erosion

    END SUBROUTINE SSC_Constant_Model
    
! - - - - - - - - - - - - - - - - 
!
!       Local Erosion
!
! - - - - - - - - - - - - - - - - 

    SUBROUTINE Local_Erosion
    
        IMPLICIT NONE
        
        INTEGER :: i , j        

#ifdef SAINT_VENENT        
        
    ! tau > critical shear stress
        
            WHERE(TAU_TD > CTAUE)

    ! Linear Function
            
                EROQ = QE0 * (TAU_TD/CTAUE-1.D0) / RHOS !* Upsilon

    ! tau <= critical shear stress
            
            ELSEWHERE(TAU_TD <= CTAUE)

    ! set erosion zero
            
                EROQ = 0.D0
            
            ENDWHERE
        
    ! The location of EROQ is not at the grid points
    ! So center average
    
            DO j = NL , NU
                DO i = ML , MU
                    IF(i==ML) THEN
                        ERO(i,j) = EROQ(i,j)
                    ELSEIF(i==MU) THEN
                        ERO(i,j) = EROQ(i-1,j)
                    ELSE
                        ERO(i,j) = 0.5D0*(EROQ(i,j)+EROQ(i-1,j))
                    END IF
                END DO
            END DO

#elif defined(CONTINUITY)

    ! tau > critical shear stress
        
            WHERE(TAU_TD > CTAUE)

    ! Linear Function
            
                ERO = QE0 * (TAU_TD/CTAUE-1.D0) / RHOS !* Upsilon

    ! tau <= critical shear stress
            
            ELSEWHERE(TAU_TD <= CTAUE)

    ! set erosion zero
            
                ERO = 0.D0
            
            ENDWHERE

#endif

    END SUBROUTINE Local_Erosion

! - - - - - - - - - - - - - - - - 
!
!       Local Total Erosion
!
! - - - - - - - - - - - - - - - - 

    SUBROUTINE Local_Total_Erosion
    
        IMPLICIT NONE
        
        INTEGER :: i
        
        DO i = ML , MU
            ERO_1D(i) = (2.D0*SUM(ERO(i,:))-ERO(i,NL)-ERO(i,NU))*dx/2.D0
        END DO
    
    END SUBROUTINE Local_Total_Erosion   
        
! - - - - - - - - - - - -         
!
!   BED LOAD TRANSPORT
!
! - - - - - - - - - - - -        
    
    SUBROUTINE  Lateral_Bedload_Transport
    
        IMPLICIT NONE
        
        INTEGER :: i , j        
        REAL(fp) :: slope1

#ifdef SAINT_VENENT
        
        DO i = ML , MU-1
            DO j = (NL+NU)/2 , NU-1
                slope1 = (DepQ(i,j+1)-DepQ(i,j))/dx
                
    ! discontinuous case
                
                !IF(TAU_TD(i,j)>CTAUB) THEN
                !    Lateral_Transport_Q(i,j) = LBFC*slope1*(TAU_TD(i,j)/CTAUB-1.D0)
                !ELSE
                !    Lateral_Transport_Q(i,j) = 0.D0
                !END IF
                
    ! continuous case
                
                Lateral_Transport_Q(i,j) = LBFC*slope1*(0.4D0*TAU_TD(i,j)**2/CTAUB/rhoRsgd)

    ! SYMMETRY
                
                Lateral_Transport_Q(i,NU-j) = - Lateral_Transport_Q(i,j) 
                
            END DO
        END DO
        
        
    ! The location of Lateral_Transport is not at the grid points
    ! So center average
    
        DO j = NL , NU-1
            DO i = ML , MU
                IF(i==ML) THEN
                    Lateral_Transport(i,j) = Lateral_Transport_Q(i,j)
                ELSEIF(i==MU) THEN
                    Lateral_Transport(i,j) = Lateral_Transport_Q(i-1,j)
                ELSE
                    Lateral_Transport(i,j) = 0.5D0*(Lateral_Transport_Q(i,j)+Lateral_Transport_Q(i-1,j))
                END IF
            END DO
        END DO 
        
#elif defined(CONTINUITY)

        DO i = ML , MU
            DO j = (NL+NU)/2 , NU-1
                slope1 = (DepQ(i,j+1)-DepQ(i,j))/dx
                
    ! discontinuous case
                
                !IF(TAU_TD(i,j)>CTAUB) THEN
                !    Lateral_Transport_Q(i,j) = LBFC*slope1*(TAU_TD(i,j)/CTAUB-1.D0)
                !ELSE
                !    Lateral_Transport_Q(i,j) = 0.D0
                !END IF
                
    ! continuous case
                
                Lateral_Transport_Q(i,j) = LBFC*slope1*(0.4D0*TAU_TD(i,j)**2/CTAUB/rhoRsgd)

    ! SYMMETRY
                
                Lateral_Transport_Q(i,NU-j) = - Lateral_Transport_Q(i,j) 
                
            END DO
        END DO    
    
#endif

    END SUBROUTINE Lateral_Bedload_Transport
    
! - - - - - - - - - - - -         
!
!   AVALANCHE
!
! - - - - - - - - - - - -  
    
    SUBROUTINE  Lateral_Avalanche_Transport
    
        IMPLICIT NONE
        
        INTEGER :: i , j        
        REAL(fp) :: slope1

! - - - - - - - - - - - -         
!
!   AVALANCHE
!
! - - - - - - - - - - - -      

        DO i = ML , MU
            DO j = NL , NU-1
                DDep(i,j) = DATAN((Dep0(i,j+1)-Dep0(i,j))/dx)
                IF(DDep(i,j)>=0.D0) THEN
                    Lateral_Avalanche(i,j) = DTAN(DMAX1(0.D0,DDep(i,j)-Repose_Angle))*Avalanche_Coef
                ELSEIF(DDep(i,j)<0.D0) THEN
                    Lateral_Avalanche(i,j) = DTAN(DMIN1(0.D0,DDep(i,j)+Repose_Angle))*Avalanche_Coef
                END IF
            END DO
        END DO        
        
    
    END SUBROUTINE Lateral_Avalanche_Transport    
    
! = = = = = = = = = = = = = = = = = = = = 
!
!   ADI ROUTINE: solve tridiagonal matrix
!
! = = = = = = = = = = = = = = = = = = = =     

    SUBROUTINE Tridiagonal_ADI(DIR,ROW)    
    
        IMPLICIT NONE
        
        CHARACTER(LEN=1) :: DIR
        INTEGER          :: ROW

        REAL(fp) :: AMT(6,3001)
        REAL(fp) :: DCE , DCW , DCN , DCS
        INTEGER  :: i , ic , k
        INTEGER  :: IDE, IDW, IDS, IDN
        ic = 0
        
        DO i = ILOW , IUP
            
            ic = ic+1

! Local point is dry

            IF(GIDP(i)==1) THEN

! set zero
                
                amt(1,ic) = 0.D0
                amt(2,ic) = 1.D0
                amt(3,ic) = 0.D0
                amt(4,ic) = 0.D0
                CYCLE
                
! Local point is wet
                
            ELSEIF(GIDP(i)==0) THEN

!- - - - - - - - - - - - - -
! Set diffusion coefficients             
!- - - - - - - - - - - - - -
                
!SOUTH
                
                IF(ROW==JLOW) THEN
                    IDS = 1
                ELSE
                    IDS = GIDS(I)
                END IF

                IF(IDS==0) THEN
                    DCS = VSCS(i)  * LT / LYY * 0.5D0 / HP(i)
                END IF  

!NORTH
                
                IF(ROW==JUP) THEN
                    IDN = 1
                ELSE
                    IDN = GIDN(I)
                END IF

                IF(IDN==0) THEN
                    DCN = VSCN(i)  * LT / LYY * 0.5D0 / HP(i)
                END IF  

!WEST

                IF(i==ILOW) THEN
                    IDW = 1
                ELSE
                    IDW = GIDP(i-1)
                END IF

                IF(IDW==0) THEN
                    DCW = VSCP(i-1) * LT / LXX * 0.5D0 / HP(i)  
                END IF  
                
!EAST                

                IF(i==IUP) THEN
                    IDE = 1
                ELSE
                    IDE = GIDP(i+1)
                END IF                
                
                IF(IDE==0) THEN
                    DCE = VSCP(i)  * LT / LXX * 0.5D0 / HP(i)
                END IF
                
!- - - - - - - - - - - - - - - - 
! Left hand side of the equation
!- - - - - - - - - - - - - - - -                 

                IF(IDE==1.AND.IDW==0) THEN
                    
                    amt(1,ic) = - 2.D0*DCW
                    amt(2,ic) = + 2.D0*DCW + 1.D0
                    amt(3,ic) = 0.D0
                    
                ELSEIF(IDW==1.AND.IDE==0) THEN
                    
                    amt(1,ic) =   0.D0
                    amt(2,ic) = + 2.D0*DCE + 1.D0
                    amt(3,ic) = - 2.D0*DCE
                    
                ELSEIF(IDW==1.AND.IDE==1) THEN
                    
                    amt(1,ic) =   0.D0
                    amt(2,ic) =   1.D0
                    amt(3,ic) =   0.D0
                    
                ELSEIF(IDW==0.AND.IDE==0)  THEN
                    
                    amt(1,ic) = - DCW
                    amt(2,ic) = + DCW + DCE + 1.D0
                    amt(3,ic) = - DCE
                    
                ELSE
                    
                    WRITE(*,*) ROW, GIDP(i-1) , GIDP(i+1)
                    PAUSE 'tri_matrix_error_1'
                    
                END IF
                
!- - - - - - - - - - - - - - - - 
! Left hand side of the equation
!- - - - - - - - - - - - - - - - 
                
                IF(    IDN==1.AND.IDS==0) THEN
                    
                    amt(4,ic) = VARP(i) + 2.D0*DCS*(VARS(i)-VARP(i))
                    
                ELSEIF(IDS==1.AND.IDN==0) THEN
                    
                    amt(4,ic) = VARP(i) + 2.D0*DCN*(VARN(i)-VARP(i))
                    
                ELSEIF(IDS==1.AND.IDN==1) THEN
                    
                    amt(4,ic) = VARP(i)
                    
                ELSEIF(IDS==0.AND.IDN==0)  THEN
                    
                    amt(4,ic) = VARP(i) + DCN*(VARN(i)-VARP(i)) + DCS*(VARS(i)-VARP(i))

                ELSE
                    
                    WRITE(*,*) ROW, IDS , IDN
                    PAUSE 'tri_matrix_error_2'
                    
                END IF
                
        ! sourse term
                
                amt(4,ic) = amt(4,ic) + CEP(i)/HP(i)*0.5D0*LT
                amt(2,ic) = amt(2,ic) + WS/HP(i)*0.5D0*LT
                
            END IF
            
        END DO

    !    OPEN(5,file='amt.txt')
    !    do i = 1, ic
    !        WRITE(5,'(10e25.10)') (amt(k,i),k=1,4)
    !    end do
    !    CLOSE(5)

        CALL solve3_LU(ic)

        fruit(:) = amt(4,:)
        
        !IF(pros=='sedi') THEN
        !IF(dir=='y') THEN
        !OPEN(5,file='amt2.txt')
        !do i = 1, ic
        !    WRITE(5,'(10e25.10)') (amt(k,i),k=1,4)
        !end do
        !CLOSE(5)
        !pause 'xg'
        !end IF

    END SUBROUTINE Tridiagonal_ADI     

    
!- - - - - - - - - - - - - -
!    
!   Bed evolution per step
!
!- - - - - - - - - - - - - -

! 2D routine
    
    SUBROUTINE Bed_Evolution(SSC1)
    
        IMPLICIT NONE
        
        REAL(fp) , DIMENSION(:,:) , POINTER :: SSC1
        INTEGER :: i , j
        
! Erosion is positive and Deposition is negative
        
        DO i = ML , MU
            DO j = NL , NU
                Evolution_Step(i,j) = ERO(i,j) - SSC1(i,j)*WS
                Evolution_TDCYC(i,j) = Evolution_TDCYC(i,j) + Evolution_Step(i,j)
            END DO
        END DO
        
    END SUBROUTINE Bed_Evolution

! 1D routine
    
    SUBROUTINE Bed_Evolution_1D(SSC1_1D,EE1)
    
        IMPLICIT NONE
        
        REAL(fp) , DIMENSION(:) , POINTER :: SSC1_1D , EE1
        REAL(fp) :: cc , MC
        REAL(fp) :: erosion0, Deposition0, bedload0, avalanche0
        INTEGER :: i , j
        

! - - - - - - - - - - - - - - - - -
!
!   SUM THE bedevolution
!   Erosion is positive and Deposition is negative
!
! - - - - - - - - - - - - - - - - - 
        
        DO i = ML , MU
            DO j = (NL+NU)/2 , NU

    ! erosion
                
                erosion0 = ERO(i,j)

    ! Deposition
                
                IF(NDD_HWL(i,j)==0) THEN
                    
#ifdef SSC_1D_DIFFUSION
                    Deposition0 = - SSC1_1D(i)*WS         
#elif defined(SSC_2D_DIFFUSION)
                    Deposition0 = - SSC1(i)*WS
#elif defined(SSC_CONSTANT)
                    Deposition0 = - SSC_CONST*WS
#endif

                ELSE
                    
                    Deposition0 = 0.D0
                    
                ENDIF

    ! bedload
                
#ifdef BEDLOAD_INC

                IF(j==NL) THEN
                    bedload0 = + Lateral_Transport(i,j)/dx
                ELSEIF(j==NU) THEN
                    bedload0 = - Lateral_Transport(i,j-1)/dx
                ELSE
                    bedload0 = (Lateral_Transport(i,j)-Lateral_Transport(i,j-1))/dx
                END IF

#else
                bedload0 = 0.D0           
#endif

    ! avalanche

#ifdef AVAlANCHE_INC

                IF(j==NL) THEN
                    avalanche0 = Lateral_Avalanche(i,j)/dx
                ELSEIF(j==NU) THEN
                    avalanche0 = - Lateral_Avalanche(i,j-1)/dx
                ELSE
                    avalanche0 = (Lateral_Avalanche(i,j)-Lateral_Avalanche(i,j-1))/dx
                END IF

#else
                avalanche0 = 0.D0           
#endif

    ! UPDATE THE BED

#ifdef CUT_OFF_SSC
                IF(CD2>CUT_OFF_SSC) THEN
                    Evolution_Step(i,j) = erosion0 + bedload0 + avalanche0
                ElSE
                    Evolution_Step(i,j) = erosion0 + Deposition0 + bedload0 + avalanche0
                ENDIF
#else
                Evolution_Step(i,j) = erosion0 + Deposition0 + bedload0 + avalanche0
#endif
                cc = Evolution_Step(i,j)*dt*MorFactor
                Dep0(i,j) = Dep0(i,j) + cc   
                Dep0(i,NL+NU-j) = Dep0(i,j)        
                
#ifndef TIDE_OSCILLATION
                IF(Dep0(i,j)<0.05D0) Dep0(i,j)=0.05D0
#endif
                
            END DO
        END DO    
        
    END SUBROUTINE Bed_Evolution_1D

! only erosion
    
    SUBROUTINE Bed_Evolution_Only_ERO(SSC1_1D)
    
        IMPLICIT NONE
        
        REAL(fp) , DIMENSION(:) , POINTER :: SSC1_1D
        INTEGER :: i , j
        
! Erosion is positive and Deposition is negative
! calculate lateral transport
        
        DO i = ML , MU
            DO j = NL , NU-1
                DDep(i,j) = DATAN((Dep1(i,j+1)-Dep1(i,j))/dx)
                IF(DDep(i,j)>=0.D0) THEN
                    Lateral_Transport(i,j) = DTAN(DMAX1(0.D0,DDep(i,j)-Repose_Angle))*Avalanche_Coef
                ELSEIF(DDep(i,j)<0.D0) THEN
                    Lateral_Transport(i,j) = DTAN(DMIN1(0.D0,DDep(i,j)+Repose_Angle))*Avalanche_Coef
                END IF
            END DO
        END DO
        
        DO i = ML , MU
            DO j = NL , NU

! erosion and Deposition
                
                Evolution_Step(i,j) = ERO(i,j)

! lateral transport
                
                IF(j==NL) THEN
                    Evolution_Step(i,j) = Evolution_Step(i,j) + Lateral_Transport(i,j)/dx
                ELSEIF(j==NU) THEN
                    Evolution_Step(i,j) = Evolution_Step(i,j) - Lateral_Transport(i,j-1)/dx
                ELSE
                    Evolution_Step(i,j) = Evolution_Step(i,j) + (Lateral_Transport(i,j)-Lateral_Transport(i,j-1))/dx
                END IF
                
                IF(CD2 > 2)  Dep1(i,j) = Dep1(i,j) + Evolution_Step(i,j)*dt*MorFactor

            END DO
        END DO    
        
    END SUBROUTINE Bed_Evolution_Only_ERO
    
! constant Deposition
    
    SUBROUTINE Bed_Evolution_Const_DepOSIT(SSC1_1D,EE1)
    
        IMPLICIT NONE
        
        REAL(fp) , DIMENSION(:) , POINTER :: SSC1_1D, EE1
        REAL(fp)                          :: c0
        INTEGER :: i , j
        
! Erosion is positive and Deposition is negative
! calculate lateral transport
        
        DO i = ML , MU
            DO j = NL , NU-1
                DDep(i,j) = DATAN((Dep0(i,j+1)-Dep0(i,j))/dx)
                IF(DDep(i,j)>=0.D0) THEN
                    Lateral_Transport(i,j) = DTAN(DMAX1(0.D0,DDep(i,j)-Repose_Angle))*Avalanche_Coef
                ELSEIF(DDep(i,j)<0.D0) THEN
                    Lateral_Transport(i,j) = DTAN(DMIN1(0.D0,DDep(i,j)+Repose_Angle))*Avalanche_Coef
                END IF
            END DO
        END DO
        
        DO i = ML , MU
            
!! sediment concentration
            
            IF(i<=30) then
                c0 = 0.D0
            ELSEIF(i>=100) then
                c0 = 0.005D0
            ELSE
                c0 = DBLE(i-30)/DBLE(100-30)*0.005D0
            END IF
            
            DO j = NL , NU

! erosion

                Evolution_Step(i,j) = ERO(i,j)

! lateral transport
                
                IF(j==NL) THEN
                    Evolution_Step(i,j) = Evolution_Step(i,j) + Lateral_Transport(i,j)/dx
                ELSEIF(j==NU) THEN
                    Evolution_Step(i,j) = Evolution_Step(i,j) - Lateral_Transport(i,j-1)/dx
                ELSE
                    Evolution_Step(i,j) = Evolution_Step(i,j) + (Lateral_Transport(i,j)-Lateral_Transport(i,j-1))/dx
                END IF

! Deposition        
        
                IF(EE1(i)+Dep0(i,j)>DDRY) THEN
                    !Evolution_Step(i,j) = Evolution_Step(i,j)*dt*MorFactor - DMIN1(0.5D0*(EE1(i)+Dep0(i,j)-DDRY),WS*c0*dt*MorFactor)
                    Evolution_Step(i,j) = Evolution_Step(i,j)*dt*MorFactor - WS*c0*dt*MorFactor
                ELSE
                    Evolution_Step(i,j) = Evolution_Step(i,j)*dt*MorFactor
                ENDIF

                IF(CD2 > 2) THEN
                    Dep0(i,j) = Dep0(i,j) + Evolution_Step(i,j)
                ENDIF
                
            END DO
        END DO    
        
    END SUBROUTINE Bed_Evolution_Const_DepOSIT
    
!= = = = = = = = 
!
! 1D sediment diffusion equation
!
!= = = = = = = =  

    SUBROUTINE Tridiagonal_1D(DIR,ROW)    
    
        IMPLICIT NONE
        
        CHARACTER(LEN=1) :: DIR
        INTEGER          :: ROW

        REAL(fp) :: DCE , DCW
        REAL(fp) :: ADE , ADW , ADP
        INTEGER  :: i , ic , k
        INTEGER  :: IDE, IDW
        ic = 0
        
        DO i = ILOW , IUP
            
            ic = ic+1

! Local point is dry

            IF(GIDP(i)==1) THEN

! set zero
                
                amt(1,ic) = 0.D0
                amt(2,ic) = 1.D0
                amt(3,ic) = 0.D0
                amt(4,ic) = 0.D0
                CYCLE
                
! Local point is wet
                
            ELSEIF(GIDP(i)==0) THEN

!- - - - - - - - - - - - - -
! Set diffusion coefficients             
!- - - - - - - - - - - - - -

!WEST

                IF(i==ILOW) THEN
                    !IDW = 1    ! NEUMANN BOUNDARY
                    IDW = -2    ! zero Dirichlet
                ELSE
                    IDW = GIDP(i-1)
                END IF

                IF(IDW==0) THEN
                    DCW = VSCP(i-1) * LT / LXX / HP(i)  
                END IF  
                
!EAST                

                IF(i==IUP) THEN
                    IDE = 1
                ELSE
                    IDE = GIDP(i+1)
                END IF                
                
                IF(IDE==0) THEN
                    DCE = VSCP(i)  * LT / LXX / HP(i)
                END IF

!- - - - - - - - - - - - - - - - 
! Residual Advection --- Upwind scheme
!- - - - - - - - - - - - - - - - 

                IF(i/=IUP.and.i/=ILOW) THEN
                    ADE = DMIN1(REQP(i),0.D0) * LT / LX 
                    ADW = -DMAX1(REQP(i-1),0.D0) * LT / LX
                    ADP = (DMAX1(REQP(i),0.D0)-DMIN1(REQP(i-1),0.D0)) * LT / LX
                ELSE
                    ADE = 0.D0
                    ADW = 0.D0
                    ADP = 0.D0
                END IF

!- - - - - - - - - - - - - - - - 
! Left hand side of the equation
!- - - - - - - - - - - - - - - -                 

                IF(IDE==1.AND.IDW==0) THEN
                    
                    amt(1,ic) = - 2.D0*DCW
                    amt(2,ic) = + 2.D0*DCW + 1.D0
                    amt(3,ic) = 0.D0
                    
                ELSEIF(IDW==1.AND.IDE==0) THEN
                    
                    amt(1,ic) =   0.D0
                    amt(2,ic) = + 2.D0*DCE + 1.D0
                    amt(3,ic) = - 2.D0*DCE
                    
                ELSEIF(IDW==1.AND.IDE==1) THEN
                    
                    amt(1,ic) =   0.D0
                    amt(2,ic) =   1.D0
                    amt(3,ic) =   0.D0
                    
                ELSEIF(IDW==0.AND.IDE==0)  THEN
                    
                    amt(1,ic) = - DCW
                    amt(2,ic) = + DCW + DCE + 1.D0
                    amt(3,ic) = - DCE
                    
                ELSEIF(IDW==-2.AND.IDE==0)  THEN   ! zero dirichlet boundary
                    
                    amt(1,ic) =   0.D0
                    amt(2,ic) = + DCW + DCE + 1.D0
                    amt(3,ic) = - DCE    
                    
                ELSE
                    
                    WRITE(*,*) ROW, GIDP(i-1) , GIDP(i+1)
                    PAUSE 'tri_matrix_error_1'
                    
                END IF

!- - - - - - - - - - - - - - - - 
! Advection term ---- RESIDUAL
!- - - - - - - - - - - - - - - - 
                
#ifdef RESIDUALINC  

                amt(1,ic) = amt(1,ic) + ADW
                amt(2,ic) = amt(2,ic) + ADP
                amt(3,ic) = amt(3,ic) + ADE

#endif

!- - - - - - - - - - - - - - - - 
! Right hand side of the equation
!- - - - - - - - - - - - - - - - 
                    
                amt(4,ic) = VARP(i)
                
        ! sourse term
                
                amt(4,ic) = amt(4,ic) + CEP(i)/HP(i)*LT
                amt(2,ic) = amt(2,ic) + WS/HP(i)*LT*BP(i)
                
            END IF
            
        END DO

    !    OPEN(5,file='amt.txt')
    !    do i = 1, ic
    !        WRITE(5,'(10e25.10)') (amt(k,i),k=1,4)
    !    end do
    !    CLOSE(5)

        CALL solve3_LU(ic)

        fruit(:) = amt(4,:)
        
        !IF(pros=='sedi') THEN
        !IF(dir=='y') THEN
        !OPEN(5,file='amt2.txt')
        !do i = 1, ic
        !    WRITE(5,'(10e25.10)') (amt(k,i),k=1,4)
        !end do
        !CLOSE(5)
        !pause 'xg'
        !end IF

    END SUBROUTINE Tridiagonal_1D    
    
    END MODULE Sediment_Transport
    
    
