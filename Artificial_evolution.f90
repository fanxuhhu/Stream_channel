#include "AMYCAI.h"

MODULE ARTIFICIAL

    USE FloatPrecision
    USE Parameters
    USE ArraySize
    USE Topography
    USE GRIDS

    CONTAINS

SUBROUTINE ARTIFICIAL_BC_MDEP 

    IMPLICIT NONE
    INTEGER :: I , J , I0 , ii
    INTEGER :: status , ISTAT
    CHARACTER(len=80) :: err_msg
    CHARACTER(len=14) :: HistoryFile
    REAL(fp) :: y, y1, y2, z1, z2, slope1, slope2, amp0
    INTEGER , DIMENSION(1) :: SEED = (/3/)
    REAL(fp) :: RANDOMDep
    INTEGER :: JJ

    BC(:) = 100.D0
    MDEP0(:) = 10.D0

END SUBROUTINE ARTIFICIAL_BC_MDEP

SUBROUTINE ARTIFICIAL_BC_MDEP_EVO(UU1)
    
    IMPLICIT NONE
    INTEGER :: I , J , I0 , ii
    INTEGER :: status , ISTAT
    CHARACTER(len=80) :: err_msg
    CHARACTER(len=14) :: HistoryFile
    REAL(fp) :: bedshear , chezy
    INTEGER , DIMENSION(1) :: SEED = (/3/)
    REAL(fp) :: RANDOMDep
    INTEGER :: JJ
    
        REAL(fp)               , DIMENSION(:) , POINTER , INTENT(IN)    :: UU1

    DO I = ML , MU-1
        chezy = 20.D0*HR(I)**(1.D0/6.D0)
        bedshear = UU1(I)**2*GRA/chezy**2*RHOW
        IF(bedshear>0.05D0) then        
            BCEVO(I) = BCEVO(I) + (bedshear-0.1D0)/0.1D0*0.001D0
        ELSE
            BCEVO(I) = BCEVO(I) - 0.5D0*0.001D0/0.05D0*bedshear
        ENDIF
    ENDDO
    BCEVO(MU) = BCEVO(MU-1)

    END SUBROUTINE ARTIFICIAL_BC_MDEP_EVO
    
END MODULE ARTIFICIAL