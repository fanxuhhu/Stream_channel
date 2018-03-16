    
SUBROUTINE Calculate_MAXDEP
    USE ArraySize
    USE Topography
    IMPLICIT NONE

    INTEGER :: i , j
    
    DO i = ML , MU
        MAXDEP(i) = MAXVAL(DEP0(i,:))
    END DO 

END SUBROUTINE Calculate_MAXDEP