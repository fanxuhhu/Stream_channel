SUBROUTINE BND_WATERLEVEL(T0,HH2)
USE FLOATPRECISION
USE PARAMETERS
USE ARRAYSIZE
USE BOUNDARY
USE DEPTH
IMPLICIT NONE
REAL(FP) , DIMENSION(:,:) , POINTER :: HH2
REAL(FP)                            :: T0
INTEGER                             :: K , I , J , I1 , I2 , J1 , J2

DO K = 1 , BNDNOS
    IF(BNDSTP(1,K)==BNDEDP(1,K)) THEN
        I = BNDSTP(1,K)
        J1 = BNDSTP(2,K)
        J2 = BNDEDP(2,K)
        DO J = J1 , J2 , 1
            HH2(I,J) = 0.5D0*SIN(T0*TDFQ) + DEP0(I,J)
        END DO
    ELSE
        J = BNDSTP(2,K)
        I1 = BNDSTP(1,K)
        I2 = BNDEDP(1,K)
        DO I = I1 , I2 , 1
            HH2(I,J) = 0.5D0*SIN(T0*TDFQ) + DEP0(I,J)
        END DO
    END IF
END DO

RETURN
END SUBROUTINE BND_WATERLEVEL