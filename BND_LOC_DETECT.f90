SUBROUTINE BND_LOC_DETECT
USE FLOATPRECISION
USE ARRAYSIZE
USE DEPTH
USE BOUNDARY
IMPLICIT NONE
INTEGER           :: M  , N
INTEGER           :: I  , J  , STATUS  , L1  , K
INTEGER           :: I1 , I2 , J1 , J2 , IJ1 , IJ2
CHARACTER(LEN=80) :: ERR_MSG

! CALCULATE THE NUMBER OF BOUNDARYS
BNDNOS = 0
L1 = 0
! ROW <<ML>>
DO J = NL , NU-1 , 1
    IF(J==NL.AND.NDD(ML,J)/=-1) THEN
        BNDNOS = BNDNOS + 1
        CYCLE
    END IF
    IF(NDD(ML,J)==-1.AND.NDD(ML,J+1)/=-1) THEN 
        BNDNOS = BNDNOS + 1
    END IF
END DO
! ROW <<MU>>
DO J = NL , NU-1 , 1
    IF(J==NL.AND.NDD(MU,J)/=-1) THEN
        BNDNOS = BNDNOS + 1
        CYCLE
    END IF
    IF(NDD(MU,J)==-1.AND.NDD(MU,J+1)/=-1) THEN 
        BNDNOS = BNDNOS + 1
    END IF
END DO
! COLUMN <<NL>>
DO I = ML , MU-1 , 1
    IF(I==ML.AND.NDD(I,NL)/=-1) THEN
        BNDNOS = BNDNOS + 1
        CYCLE
    END IF
    IF(NDD(I,NL)==-1.AND.NDD(I+1,NL)/=-1) THEN 
        BNDNOS = BNDNOS + 1
    END IF
END DO  
! COLUMN <<NU>>  
DO I = ML , MU-1 , 1
    IF(I==ML.AND.NDD(I,NU)/=-1) THEN
        BNDNOS = BNDNOS + 1
        CYCLE
    END IF
    IF(NDD(I,NU)==-1.AND.NDD(I+1,NU)/=-1) THEN 
        BNDNOS = BNDNOS + 1
    END IF
END DO

! DETECT THE END POINTS OF EACH BOUNDARY
ALLOCATE(BNDSTP(2,BNDNOS),STAT=STATUS,ERRMSG=ERR_MSG)
ALLOCATE(BNDEDP(2,BNDNOS),STAT=STATUS,ERRMSG=ERR_MSG)
IJ1 = 0
IJ2 = 0
! ROW <<ML>>
M = ML
N = 0
CALL GET_BNDSTP(NDD,IJ1,ML,MU,NL,NU,M,N)
CALL GET_BNDEDP(NDD,IJ2,ML,MU,NL,NU,M,N)
! ROW <<MU>>
M = MU
N = 0
CALL GET_BNDSTP(NDD,IJ1,ML,MU,NL,NU,M,N)
CALL GET_BNDEDP(NDD,IJ2,ML,MU,NL,NU,M,N)
! COLUMN <<NL>>
M = 0
N = NL
CALL GET_BNDSTP(NDD,IJ1,ML,MU,NL,NU,M,N)
CALL GET_BNDEDP(NDD,IJ2,ML,MU,NL,NU,M,N)
! COLUMN <<NU>>
M = 0
N = NU
CALL GET_BNDSTP(NDD,IJ1,ML,MU,NL,NU,M,N)
CALL GET_BNDEDP(NDD,IJ2,ML,MU,NL,NU,M,N)

RETURN
END SUBROUTINE BND_LOC_DETECT


SUBROUTINE GET_BNDSTP(NDD,IJ1,ML,MU,NL,NU,M,N)
USE BOUNDARY
IMPLICIT NONE
INTEGER :: ML , MU , NL , NU , M , N
INTEGER , DIMENSION(ML:MU,NL:NU) :: NDD
INTEGER :: IJ1
INTEGER :: I , J 

IF(M/=0.AND.N==0) THEN
    DO J = NL , NU-1 , 1
        IF(J==NL.AND.NDD(M,J)/=-1) THEN
            IJ1 = IJ1 + 1
            BNDSTP(1,IJ1) = M
            BNDSTP(2,IJ1) = J
            CYCLE
        END IF
        IF(NDD(M,J)==-1.AND.NDD(M,J+1)/=-1) THEN 
            IJ1 = IJ1 + 1
            BNDSTP(1,IJ1) = M
            BNDSTP(2,IJ1) = J+1
            CYCLE
        END IF
    END DO
ELSE IF(M==0.AND.N/=0) THEN
    DO I = ML , MU-1 , 1
        IF(I==ML.AND.NDD(I,N)/=-1) THEN
            IJ1 = IJ1 + 1
            BNDSTP(1,IJ1) = I
            BNDSTP(2,IJ1) = N
            CYCLE
        END IF
        IF(NDD(I,N)==-1.AND.NDD(I+1,N)/=-1) THEN 
            IJ1 = IJ1 + 1
            BNDSTP(1,IJ1) = I+1
            BNDSTP(2,IJ1) = N
            CYCLE
        END IF
    END DO  
END IF

RETURN
END SUBROUTINE GET_BNDSTP

SUBROUTINE GET_BNDEDP(NDD,IJ2,ML,MU,NL,NU,M,N)
USE BOUNDARY
IMPLICIT NONE
INTEGER :: ML , MU , NL , NU , M , N
INTEGER , DIMENSION(ML:MU,NL:NU) :: NDD
INTEGER :: IJ2
INTEGER :: I , J 

IF(M/=0.AND.N==0) THEN
    DO J = NL+1 , NU , 1
        IF(J==NU.AND.NDD(M,J)/=-1) THEN
            IJ2 = IJ2 + 1
            BNDEDP(1,IJ2) = M
            BNDEDP(2,IJ2) = J
            CYCLE
        END IF
        IF(NDD(M,J)==-1.AND.NDD(M,J-1)/=-1) THEN 
            IJ2 = IJ2 + 1
            BNDEDP(1,IJ2) = M
            BNDEDP(2,IJ2) = J-1
            CYCLE
        END IF
    END DO
ELSE IF(M==0.AND.N/=0) THEN
    DO I = ML+1 , MU , 1
        IF(I==MU.AND.NDD(I,N)/=-1) THEN
            IJ2 = IJ2 + 1
            BNDEDP(1,IJ2) = I
            BNDEDP(2,IJ2) = N
            CYCLE
        END IF
        IF(NDD(I,N)==-1.AND.NDD(I-1,N)/=-1) THEN 
            IJ2 = IJ2 + 1
            BNDEDP(1,IJ2) = I-1
            BNDEDP(2,IJ2) = N
            CYCLE
        END IF
    END DO  
END IF

RETURN
END SUBROUTINE GET_BNDEDP