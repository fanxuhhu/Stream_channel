#include "AMYCAI.h"

MODULE Geology_Information

    USE TimeCatcher
    USE FloatPrecision
    USE ArraySize
    USE GRIDS
    USE Parameters

    IMPLICIT NONE
    
    CONTAINS

! = = = = = = = = = = = = = = = = =
!
! assignment to Seabed Parameters
!
! = = = = = = = = = = = = = = = = =
    
    SUBROUTINE BedParameters( DEP0, MDEP0, DDEP, Upsilon, MAREA)

        IMPLICIT NONE
        
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)   ,          INTENT(IN)    :: DEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)     ,          INTENT(INOUT) :: MDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)   ,          INTENT(INOUT) :: DDEP
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)   , TARGET , INTENT(INOUT) :: Upsilon
        REAL(fp) , ALLOCATABLE , DIMENSION(:)     ,          INTENT(INOUT) :: MAREA
        
        INTEGER           :: I , J , K
        INTEGER           :: status
        CHARACTER(len=80) :: err_msg

! - - - - - - - - - - - - - - - - - -
!
!   Calculate the mean water Depth
!
! - - - - - - - - - - - - - - - - - -
        
        MDEP0(:) = (2.D0*SUM(DEP0,DIM=2) - DEP0(:,NL) - DEP0(:,NU)) / 2.D0 * DX / BT

! - - - - - - - - - - - - - - - - - -
!
!   Calculate the MeanArea
!
! - - - - - - - - - - - - - - - - - -
        
        DO I = ML , MU
            MAREA(I) = 0.D0
            DO J = NL , NU-1
                MAREA(I) = MAREA(I) + (DMAX1(DEP0(I,J),DDRY)+DMAX1(DEP0(I,J+1),DDRY))/2.D0*DX
            ENDDO
        ENDDO
            
        
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!   Get the lateral slope and the factors of the shear stress (Upsilon)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
        DO I = ML , MU
            DO J = NL , NU
                IF(J==NU) THEN
                    DDEP(I,J) = 0.D0
                ELSE
                    DDEP(I,J) = (DEP0(I,J+1)-DEP0(I,J)) / dx
                END IF
                Upsilon(I,J) = SQRT(1+DDEP(I,J)**2)
            END DO     
        END DO

! sediment stratification
                
        
    END SUBROUTINE BedParameters

! = = = = = = = = = = = = = = = = =
!
! assignment to Sedimentation Parameters
! These parameters serve the sediment transport
! Model
!
! = = = = = = = = = = = = = = = = =            
            
    SUBROUTINE Sedimentation_Parameters(DEP0,TDEP0,HTDEP0,MDEP0,HMTDEP0,HETA,HND,HNDD,HDIFF)

        IMPLICIT NONE

        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,          INTENT(IN)    :: DEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,          INTENT(IN)    :: TDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,          INTENT(INOUT) :: HTDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,          INTENT(IN)    :: MDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,          INTENT(INOUT) :: HMTDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,          INTENT(INOUT) :: HETA
        INTEGER  , ALLOCATABLE , DIMENSION(:)   ,          INTENT(INOUT) :: HND
        INTEGER  , ALLOCATABLE , DIMENSION(:,:) ,          INTENT(INOUT) :: HNDD
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,          INTENT(INOUT) :: HDIFF
        !---------------------------------------------------------------------------
        INTEGER                                                          :: I , J
        INTEGER                                                          :: status
        CHARACTER(LEN=80)                                                :: err_msg 
    
#ifndef SSC_CONSTANT

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! get the water Depth at high tide
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -    

        DO I = ML , MU
            HTDEP0(I,:) = DEP0(I,:) + TDAMP0
        END DO

    ! - - - - - - - - - - - - - 
    !   get the HWL Mean Total water depth
    ! - - - - - - - - - - - - - 
    
        HMTDEP0(:) = DEP0(:,NU) + TDAMP0
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! get the wet AND dry at high tide
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
        
        WHERE(HTDEP0 > 0.D0)        ! wet  
            HNDD = 0            
        ELSEWHERE(HTDEP0 <= 0.D0)   ! dry
            HNDD = 1        
        END WHERE

        WHERE(HMTDEP0 > 0.D0)        ! wet  
            HND = 0            
        ELSEWHERE(HMTDEP0 <= 0.D0)   ! dry
            HND = 1        
        END WHERE    
        HND(ML) = -5

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  Get the diffusion coefficient for SSC transport.
    !  The locations of the diffusion coefficients are between
    !  grid points... 1D
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!    ! diffusion in y direcion (Longitudinal)
!    
!        DO I = ML , MU-1
!            IF(HMTDEP0(I)>0.D0.AND.HMTDEP0(I+1)>0.D0) THEN
!                HDIFF(I) = TD_Excursion_Y*0.5D0*(HMTDEP0(I)+HMTDEP0(I+1))
!            ELSE
!                HDIFF(I) = 0.D0
!            END IF
!        END DO
!    
!    ! Reset Maximum eta
!
!        IF(T0>=TDCYC) THEN
!            HETA(:)  = 0.D0
!        END IF
!        
#else

        WHERE(TDEP0 > 0.D0)        ! wet  
            HNDD = 0            
        ELSEWHERE(TDEP0 <= 0.D0)   ! dry
            HNDD = 1        
        END WHERE

#endif

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  Get the diffusion coefficient for SSC transport.
    !  The locations of the diffusion coefficients are between
    !  grid points... 2D
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -       

    ! diffusion in x direcion (transversal)

        !DO I = ML , MU
        !    DO J = NL , NU-1
        !        IF(HTDEP0(I,J)>0.D0.AND.HTDEP0(I,J+1)>0.D0) THEN
        !            SSC_DIFFU_HWL_X(I,J) = TD_Excursion_X*0.5D0*(HTDEP0(I,J)+HTDEP0(I,J+1))
        !        ELSE
        !            SSC_DIFFU_HWL_X(I,J) = 0.D0
        !        END IF
        !    END DO
        !END DO    

    ! diffusion in y direcion (Longitudinal)
    
        !DO I = ML , MU-1
        !    DO J = NL , NU
        !        IF(HTDEP0(I,J)>0.D0.AND.HTDEP0(I+1,J)>0.D0) THEN
        !            SSC_DIFFU_HWL_Y(I,J) = TD_Excursion_Y*0.5D0*(HTDEP0(I,J)+HTDEP0(I+1,J))
        !        ELSE
        !            SSC_DIFFU_HWL_Y(I,J) = 0.D0
        !        END IF
        !    END DO
        !END DO
    
    END SUBROUTINE Sedimentation_Parameters
            
END MODULE Geology_Information