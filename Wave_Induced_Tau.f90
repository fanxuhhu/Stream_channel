#include "AMYCAI.h"

MODULE Wave_Effects

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
!   Wave induced bed shear stress
!
!= = = = = = = = = = = = = = = = = = = = = =  
    
    !SUBROUTINE Wave_Induced_Tau(WVTAU,TDEP0)
    !
    !    IMPLICIT NONE
    !    
    !    REAL(fp)                                          , INTENT(IN)    :: T0       
    !    REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: UU1
    !    REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: UU2
    !    REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: EE1 
    !    REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: EE2   
    !    REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: QQ1
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(IN)    :: BC
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: HR
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET  , INTENT(IN)    :: FC
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET  , INTENT(INOUT) :: HETA
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET  , INTENT(INOUT) :: Q_RESIDUAL2
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: DEP0
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(INOUT) :: TDEP0
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(IN)    :: MDEP0
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(IN)    :: MAREA
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: MDEPE
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: MTDEPE
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: EtaU
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: MDEPU
    !    REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: MTDEPU
    !    INTEGER  , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: ND
    !    INTEGER  , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: NDQ
    !    
    !    REAL(fp) ,               DIMENSION(:) , POINTER                 :: TEMP
    !    LOGICAL                                                         :: GSTB
    !    INTEGER                                                         :: I , J
    !    REAL(fp)                                                        :: A1,A2,B1,B2
    !        
    !    DO I = ML , MU
    !        DO J = NL , NU
    !        A1 = 0.493D0*delta**(0.75D0)
    !        A2 = 0.331D0*delta**(1.01D0)
    !        B1 = 0.00313*kai**(0.57D0)
    !        B2 = 5.215D-4*kai**(0.73)
    !        epsilon = 0.00364D0*(TANH(A1)*TANH(B1/TANH(A1)))**(1.74D0)
    !        zeta = 0.133D0*(TANH(A2)*TANH(B2/TANH(A2)))**(-0.37D0)
    !        ENDDO
    !    ENDDO
    !        
    !        
    !END SUBROUTINE Wave_Induced_Tau

    
! - - - - - - - -

END MODULE Wave_Effects


    
    