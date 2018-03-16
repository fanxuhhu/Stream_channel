! This routine is used to decide whether the cross-sections are
! wet OR dry, AND give them corrsponding ID:
!  0 <=> wet
!  1 <=> dry
! -2 <=> OPEN boundary
!  2 <=> closed boundary 
    
#include "AMYCAI.h"

SUBROUTINE GRID_ND

    USE FloatPrecision
    USE ArraySize
    USE hydro
    USE Topography
    USE Parameters

    IMPLICIT NONE
    
    REAL(fp) :: MAXH , DDRY2
    INTEGER  :: I

    NDQ(ML)  = -2 ! OPEN boundary   
    NDQ(MU) =  2 ! CLOSE boundary
    
    DO I = ML+1 , MU-1
        IF(MTDEPU(I) > DDRY) THEN
            NDQ(I) = 0
        ELSE
            NDQ(I) = 1
        END IF
    END DO
    
    ND(:) = 0
    ND(ML) = -2
    ND(MU+1) = 2
    DO I = ML+1 , MU
        IF(MTDEPE(I) < DDRY) THEN
            ETA1(I) = DDRY - MDEPE(I)
        ENDIF
    END DO
    ETA1(ML) = ETA1(ML+1)
    ETA1(MU+1) = ETA1(MU)
    
END SUBROUTINE GRID_ND