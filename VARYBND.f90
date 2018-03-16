SUBROUTINE VARYBND(SURFACE,MAXDEP,DISCHARGE,ND,NDQ,GSTB)
    
USE TimeCatcher
USE FloatPrecision
USE ArraySize
USE Parameters

IMPLICIT NONE

REAL(fp) , DIMENSION(ML:MU)       :: SURFACE , DISCHARGE
REAL(fp) , DIMENSION(ML:MU)       :: MAXDEP
INTEGER  , DIMENSION(ML:MU)       :: ND , NDQ
LOGICAL  :: GSTB

REAL(fp) :: MAXH , AVH
INTEGER  :: I
INTEGER  :: CNODES


    GSTB = .TRUE.
    CNODES = 0

!   IF there are negative water depth (1D)
    
    DO I = ML , MU , 1
        MAXH   = MAXDEP(I) + SURFACE(I)
        IF(MAXH<=0.D0) THEN
            ND(I) = 1
            GSTB = .FALSE.
            CNODES = CNODES+1
        END IF
    END DO
    if(cd==720) write(*,*) CNODES
    IF(.NOT.GSTB) RETURN

    DO I = ML , MU , 1
        IF(ND(I)==1) THEN
            ND(I) = 0
        END IF
    END DO


! CHECK NDQ
    
    DO I = ML , MU-1 , 1
        IF(DISCHARGE(I)>0.D0) THEN
            AVH = MAXDEP(I) + SURFACE(I)
        ELSEIF(DISCHARGE(I)<0.D0) THEN
            AVH = MAXDEP(I+1) + SURFACE(I+1)
        ELSEIF(DISCHARGE(I)==0.D0) THEN
            AVH = 0.5D0*(MAXDEP(I)+SURFACE(I)+MAXDEP(I+1)+SURFACE(I+1))
        END IF
        IF(NDQ(I)==0) THEN
            IF(AVH<=DDRY)  THEN
                NDQ(I) = 1
            END IF
        ELSEIF(NDQ(I)==1) THEN
            IF(AVH>=DWET) NDQ(I) = 0
        END IF
    END DO

END SUBROUTINE VARYBND

