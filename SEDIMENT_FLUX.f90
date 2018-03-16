SUBROUTINE SEDI_FLUX(DIR,HH2,PP1,PP2,SXP,QQ1,SXQ)
    USE TEST
    USE FLOATPRECISION
    USE PARAMETERS
    USE DEPTH
    USE ARRAYSIZE
    IMPLICIT NONE
    CHARACTER*1                         :: DIR
    REAL(FP) , DIMENSION(:,:) , POINTER :: HH2 , PP1 , PP2 , SXP , QQ1 , SXQ
    INTEGER                             :: I , J , I1 , I2 , J1 , J2
    REAL(FP)                            :: H0 , H1, P0 , Q0 , POW1 , MAN3
    
    POW1 = 11.D0/2.D0
    MAN3 = MANNING**3
    
    IF(DIR=='X') THEN
        
        DO I = ML+1 , MU-1
            LPJ1:DO J = NL , NU-1
                P0 = 0.5D0*(PP1(I,J)+PP2(I,J))
                IF(    P0>0.D0) THEN
                    H0 = HH2(I,J)
                    H1 = HH2(I,J+1)
                ELSEIF(P0<0.D0) THEN
                    H0 = HH2(I,J+1)
                    H1 = HH2(I,J)
                ELSEIF(P0==0.D0)THEN
                    SXP(I,J) = 0.D0
                    CYCLE LPJ1
                END IF
                IF(H0>=DCRI.AND.H1>=DCRI) THEN
                    SXP(I,J) = 0.05D0*P0**5/SQRT(GRA)/MAN3/H0**POW1/REDY**2/D50
                ELSEIF(H0<DCRI.OR.H1<DCRI) THEN
                    SXP(I,J) = 0.D0
                ELSE
                    WRITE(*,*) 'ERROR FROM SEDIMENT_FLUX SUBROUTINE'
                    PAUSE
                END IF
                IF(DABS(SXP(I,J))>100.D0) THEN
                    WRITE(*,*) 'p' , I , J
                    WRITE(*,*) CD2
                    WRITE(*,*) P0 , H0
                    WRITE(*,*) CD2
                    WRITE(*,*) DIR
                    EVOGO = .FALSE.
                    CALL OUTPUT('CRI_','.DEP',CD2,DEP0,ML,MU,1,NL,NU,1)
                    PAUSE
                END IF
            END DO LPJ1
        END DO
        ! 136906
        DO J = NL+1 , NU-1
            LPJ2:DO I = ML , MU-1
                Q0 = QQ1(I,J)
                IF(    Q0>0.D0) THEN
                    H0 = HH2(I,J)
                    H1 = HH2(I+1,J)
                ELSEIF(Q0<0.D0) THEN
                    H0 = HH2(I+1,J)
                    H1 = HH2(I,J)
                ELSEIF(Q0==0.D0)THEN
                    SXQ(I,J) = 0.D0
                    CYCLE LPJ2
                END IF
                IF(H0>=DCRI.AND.H1>=DCRI) THEN
                    SXQ(I,J) = 0.05D0*Q0**5/SQRT(GRA)/MAN3/H0**POW1/REDY**2/D50
                ELSEIF(H0<DCRI.OR.H1<DCRI) THEN
                    SXQ(I,J) = 0.D0
                ELSE
                    WRITE(*,*) 'ERROR FROM SEDIMENT_FLUX SUBROUTINE'
                    PAUSE
                END IF    
                IF(DABS(SXQ(I,J))>100.D0) THEN
                    WRITE(*,*) 'q',I , J
                    WRITE(*,*) CD2
                    WRITE(*,*) Q0 , H0
                    WRITE(*,*) CD2
                    WRITE(*,*) DIR
                    EVOGO = .FALSE.
                    CALL OUTPUT('CRI_','.DEP',CD2,DEP0,ML,MU,1,NL,NU,1)
                    PAUSE
                END IF
            END DO LPJ2
        END DO
        
    ELSEIF(DIR=='Y') THEN
        
        DO I = ML+1 , MU-1
            LPJ3:DO J = NL , NU-1
                Q0 = QQ1(I,J)
                IF(    Q0>0.D0) THEN
                    H0 = HH2(I,J)
                    H1 = HH2(I,J+1)
                ELSEIF(Q0<0.D0) THEN
                    H0 = HH2(I,J+1)
                    H1 = HH2(I,J)
                ELSEIF(Q0==0.D0)THEN
                    SXQ(I,J) = 0.D0
                    CYCLE LPJ3
                END IF
                IF(H0>=DCRI.AND.H1>=DCRI) THEN
                    SXQ(I,J) = 0.05D0*Q0**5/SQRT(GRA)/MAN3/H0**POW1/REDY**2/D50
                ELSEIF(H0<DCRI.OR.H1<DCRI) THEN
                    SXQ(I,J) = 0.D0
                ELSE
                    WRITE(*,*) 'ERROR FROM SEDIMENT_FLUX SUBROUTINE'
                    PAUSE
                END IF
                IF(DABS(SXQ(I,J))>100.D0) THEN
                    WRITE(*,*) 'p' , I , J
                    WRITE(*,*) CD2
                    WRITE(*,*) Q0 , H0
                    WRITE(*,*) CD2
                    WRITE(*,*) DIR
                    EVOGO = .FALSE.
                    CALL OUTPUT('CRI_','.DEP',CD2,DEP0,ML,MU,1,NL,NU,1)
                    PAUSE
                END IF
            END DO LPJ3
        END DO

        DO J = NL+1 , NU-1
            LPJ4:DO I = ML , MU-1
                P0 = 0.5D0*(PP2(I,J)+PP1(I,J))
                IF(    P0>0.D0) THEN
                    H0 = HH2(I,J)
                    H1 = HH2(I+1,J)
                ELSEIF(P0<0.D0) THEN
                    H0 = HH2(I+1,J)
                    H1 = HH2(I,J)
                ELSEIF(P0==0.D0)THEN
                    SXP(I,J) = 0.D0
                    CYCLE LPJ4
                END IF
                IF(H0>=DCRI.AND.H1>=DCRI) THEN
                    SXP(I,J) = 0.05D0*P0**5/SQRT(GRA)/MAN3/H0**POW1/REDY**2/D50
                ELSEIF(H0<DCRI.OR.H1<DCRI) THEN
                    SXP(I,J) = 0.D0
                ELSE
                    WRITE(*,*) 'ERROR FROM SEDIMENT_FLUX SUBROUTINE'
                    PAUSE
                END IF
               IF(DABS(SXP(I,J))>100.D0) THEN
                    WRITE(*,*) 'q', I , J
                    WRITE(*,*) CD2
                    WRITE(*,*) P0 , H0
                    WRITE(*,*) CD2
                    WRITE(*,*) DIR
                    EVOGO = .FALSE.
                    CALL OUTPUT('CRI_','.DEP',CD2,DEP0,ML,MU,1,NL,NU,1)
                    PAUSE
                END IF
            END DO LPJ4
        END DO
        
    END IF
    
END SUBROUTINE SEDI_FLUX