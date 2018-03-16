SUBROUTINE Bed_Evolution(SSC1)
    USE TEST
    USE FLOATPRECISION
    USE PARAMETERS
    USE DEPTH
    USE GRIDS
    USE ARRAYSIZE    
    IMPLICIT NONE
    REAL(FP) , DIMENSION(:,:) , POINTER :: SXP , SXQ
    REAL(FP) :: DT05
    INTEGER :: I , J
    INTEGER :: I0 , J0
    DT05 = 0.5D0*DT
    
    DO I = ML+1 , MU-1 , 1
        DO J = NL+1 , NU-1 , 1
            IF(NDD(I,J)<=0.D0) THEN
                EVO(I,J) = + (SXP(I,J)-SXP(I,J-1))/DX*DT05 &
                         & + (SXQ(I,J)-SXQ(I-1,J))/DY*DT05
                IF(I>30) THEN
                    DEP0(I,J) = DEP0(I,J) + EVO(I,J)*MORACC
                END IF
                IF(DABS(EVO(I,J))>MAXEVO) THEN
                    MAXEVO = EVO(I,J)
                    MEI = I
                    MEJ = J
                END IF
            END IF
        END DO
    END DO
    
END SUBROUTINE BED_EVOLUTION