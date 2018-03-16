SUBROUTINE GRID_IDN_JDN
USE ARRAYSIZE
USE DEPTH
IMPLICIT NONE
INTEGER :: I , J , ISTAT

DO I = ML , MU , 1
    LPJ1:DO J = NL , NU-1 , 1
        IF(NDD(I,J)/=-1) THEN
            IDN(1,I) = J
            EXIT LPJ1
        END IF
    END DO LPJ1
    LPJ2:DO J = NU-1 , NL , -1
        IF(NDD(I,J)/=-1) THEN
            IDN(2,I) = J
            EXIT LPJ2
        END IF
    END DO LPJ2
END DO

DO J = NL , NU , 1
    LPJ3:DO I = ML , MU-1 , 1
        IF(NDD(I,J)/=-1) THEN
            JDN(1,J) = I
            EXIT LPJ3
        END IF
    END DO LPJ3
    LPJ4:DO I = MU-1 , ML , -1
        IF(NDD(I,J)/=-1) THEN
            JDN(2,J) = I
            EXIT LPJ4
        END IF
    END DO LPJ4
END DO

OPEN(UNIT=10,FILE='IDN.TXT',STATUS='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
DO I = ML , MU , 1
    WRITE(10,'(3I8)') I,(IDN(J,I),J=1,2)
END DO
CLOSE(10)
OPEN(UNIT=11,FILE='JDN.TXT',STATUS='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
DO J = NL , NU , 1
    WRITE(11,'(3I8)') J,(JDN(I,J),I=1,2)
END DO
CLOSE(11)

END SUBROUTINE GRID_IDN_JDN