SUBROUTINE GRID_NDD_PROGDET
    USE ARRAYSIZE
    USE DEPTH
    IMPLICIT NONE
    INTEGER :: I , J , ISTAT

    NDD(:,:)  = 0
    NDD(ML,:) = -2
!    NDD(MU,:) = -2
!    NDD(:,NL) = -2
!    NDD(:,NU) = -2
    DNDD(:,:) = 0
    DO I = ML , MU , 1
        DO J = NL , NU , 1
            IF(DEP0(I,J)==-99.99D0) NDD(I,J) = 2
        END DO
    END DO

    OPEN(UNIT=10,FILE='NDD.TXT' ,STATUS='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
    DO I = MU , ML , -1
        WRITE(10,'(500I3)')   (NDD(I,J) ,J=NL,NU)
    END DO
    CLOSE(10)
    CLOSE(11)

    END SUBROUTINE GRID_NDD_PROGDET


    SUBROUTINE GRID_NDD_USERSET
    USE ARRAYSIZE
    USE DEPTH
    IMPLICIT NONE
    INTEGER :: I , J , ISTAT

    NDD(:,:) = 0
    NDD(ML,:) = -2
    NDD(MU,:) = -2
    NDD(:,NL) = -2
    NDD(:,NU) = -2
    DO I = ML , MU , 1
        DO J = NL , NU , 1
            IF(DEP0(I,J)==-99.99D0) NDD(I,J) = -1
        END DO
    END DO

END SUBROUTINE GRID_NDD_USERSET