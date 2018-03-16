SUBROUTINE CONTINUITY_MOMENTUM_Y(PP1,PP2,HH1,HH2,QQ1,QQ2,T0,DT,DX,DY,I1,I2,IP,UPS)
USE TEST
USE FLOATPRECISION
USE PARAMETERS
USE DEPTH
IMPLICIT NONE
REAL(FP) , DIMENSION(:,:) , POINTER :: PP1 , PP2 , HH1 , HH2 , QQ1 , QQ2
REAL(FP)                            :: T0  , DT  , DX  , DY
INTEGER                             :: I1  , I2  , IP
LOGICAL                             :: UPS

REAL(FP) :: FROUDE , WF
REAL(FP) :: ADVE , ADVW
REAL(FP) :: VN , VS , PN , PS , CADL , CADR , CDIFFL , CDIFFR
REAL(FP) :: RESIST , CHEZY , QAV , HAV , GRV
REAL(FP) :: DIFFX , DIFFY
REAL(FP) :: HW , HE , HN , HS , PHIN , PHIS
INTEGER  :: ML , MU , NL , NU
INTEGER  :: I , J , K , JC , JJC , J0
REAL(FP) :: AMT(6,3001)
REAL(FP) :: SEABND
REAL(FP) :: DTDY , DTDX , DXDY , DTDXDY , DTDXDDY , DTDYDDX , DXDY4 , DT0125 , DT2 , DTDY2
REAL(FP) :: POW1 , GRA05 , GRA0125 , GDM2 , COE1
CHARACTER(LEN=80) :: ERR_MSG
INTEGER :: STATUS

DT2     = 2.D0*DT
DT0125  = 0.125D0*DT
DTDX    = DT*DX
DTDY    = DT*DY
DTDY2   = 2.D0*DTDY
DXDY    = DX*DY
DXDY4   = 4.D0*DXDY
DTDXDY  = DT*DX*DY
DTDXDDY = DT*DX/DY
DTDYDDX = DT*DY/DX
GRA05   = 0.5D0*GRA
GRA0125 = 0.125D0*GRA
GDM2    = GRA/MANNING**2
COE1    = 4.D0/3.D0
POW1    = 7.D0/3.D0

ML = LBOUND(PP1,1) 
MU = UBOUND(PP1,1)
NL = LBOUND(PP1,2) 
NU = UBOUND(PP1,2)

LPI0:DO I = I1 , I2 , IP
    ! CONTINUITY EQUATION
    JC = -1
    LPJ0:DO J = ML , MU , 1
        JC = JC + 2
        IF(NDD(J,I)>0) THEN
            AMT(1,JC) = 0.D0
            AMT(2,JC) = 0.D0
            AMT(3,JC) = 1.D0
            AMT(4,JC) = 0.D0
            AMT(5,JC) = 0.D0
            AMT(6,JC) = HH1(J,I) 
            CYCLE LPJ0   
        END IF
        IF(NDD(J,I)==-2) THEN
            AMT(1,JC) = 0.D0
            AMT(2,JC) = 0.D0
            AMT(3,JC) = 1.D0
            AMT(4,JC) = 0.D0
            AMT(5,JC) = 0.D0
            AMT(6,JC) = HH2(J,I)
            CYCLE LPJ0
        END IF
        AMT(1,JC) = 0.D0
        AMT(2,JC) = - DTDY
        AMT(3,JC) = + DXDY4                                   
        AMT(4,JC) = + DTDY
        AMT(5,JC) = 0.D0
        AMT(6,JC) = + DXDY4*HH1(J,I) &                       
                  & - DTDX*(QQ2(J,I)-QQ2(J,I-1)+QQ1(J,I)-QQ1(J,I-1)) &
                  & - DTDY*(PP1(J,I)-PP1(J-1,I))
    END DO LPJ0

    ! MOMENTUM EQUATION
    JJC = 0
    LPJ1:DO J = ML , MU-1 , 1
        JJC = JJC + 2 ! ID OF BAND MATRIX

        ! CLOSE BOUNDARY CONDITION
        IF(NDD(J,I)>0.OR.NDD(J+1,I)>0.OR.NDDQ(J,I)>0) THEN
            AMT(1,JJC) = 0.D0
            AMT(2,JJC) = 0.D0
            AMT(3,JJC) = 1.D0
            AMT(4,JJC) = 0.D0
            AMT(5,JJC) = 0.D0
            AMT(6,JJC) = 0.D0
            CYCLE LPJ1
        END IF

        !  OPEN BOUNDARY CONDITION --- SURFACE LEVEL
        IF(NDD(J,I)==-2.OR.NDD(J+1,I)==-2) THEN
                GRV  = GRA05*(HH1(J,I)+HH1(J+1,I))
                AMT(1,JJC) = 0.D0
                AMT(2,JJC) = - GRV*DTDY
                AMT(3,JJC) = + DXDY
                AMT(4,JJC) = + GRV*DTDY
                AMT(5,JJC) = 0.D0
                AMT(6,JJC) = + DXDY*PP1(J,I) &
                           & + GRV*DTDY*(DEP0(J+1,I)-DEP0(J,I))  
            CYCLE LPJ1                   
        END IF
        
        ! CROSS MOMENTUM AND CROSS DIFFUSION
        IF(NDD(J,I+1)==0.AND.NDD(J+1,I+1)==0) THEN
            HN =0.25D0*(HH1(J,I)+HH1(J+1,I)+HH1(J,I+1)+HH1(J+1,I+1))
            PN = PP1(J,I+1)
        ELSE IF(NDD(J,I+1)/=0.AND.NDD(J+1,I+1)==0) THEN
            HN = 0.25D0*(HH1(J,I)+HH1(J+1,I)+HH1(J,I)+HH1(J+1,I+1))
            PN = PP1(J,I)
        ELSE IF(NDD(J,I+1)==0.AND.NDD(J+1,I+1)/=0) THEN
            HN = 0.25D0*(HH1(J,I)+HH1(J+1,I)+HH1(J,I+1)+HH1(J+1,I))
            PN = PP1(J,I)
        ELSE
            HN = 0.5D0*(HH1(J,I)+HH1(J+1,I))
            PN = PP1(J,I)
        END IF
        IF(NDD(J,I-1)==0.AND.NDD(J+1,I-1)==0) THEN
            HS = 0.25D0*(HH1(J,I)+HH1(J+1,I)+HH1(J,I-1)+HH1(J+1,I-1))
            PS = PP1(J,I-1)
        ELSE IF(NDD(J,I-1)/=0.AND.NDD(J+1,I-1)==0) THEN
            HS = 0.25D0*(HH1(J,I)+HH1(J+1,I)+HH1(J,I)+HH1(J+1,I-1))
            PS = PP1(J,I)
        ELSE IF(NDD(J,I-1)==0.AND.NDD(J+1,I-1)/=0) THEN
            HS = 0.25D0*(HH1(J,I)+HH1(J+1,I)+HH1(J,I-1)+HH1(J+1,I))
            PS = PP1(J,I)
        ELSE
            HS = 0.5D0*(HH1(J,I)+HH1(J+1,I))
            PS = PP1(J,I)
        END IF

        ! GRAVITY 
        GRV  = GRA05*(HH1(J+1,I)+HH1(J,I))
        
        ! RESISTANCE TERM
        QAV = 0.125D0*(QQ1(J,I-1)+QQ1(J+1,I-1)+QQ2(J,I-1)+QQ2(J+1,I-1)+QQ1(J,I)+QQ1(J+1,I)+QQ2(J,I)+QQ2(J+1,I))
        IF(PP1(J,I)>=0.D0) THEN
            HAV = HH1(J,I)
        ELSE IF(PP1(J,I)<0.D0) THEN
            HAV = HH1(J+1,I)
        END IF
        RESIST = GDM2*SQRT(PP1(J,I)**2+QAV**2)/HAV**POW1

        PHIN = MANNING*SQRT(GRA)*HN**(-1.D0/6.D0)*0.4D0*0.5D0*DABS(PN+PP1(J,I))
        PHIS = MANNING*SQRT(GRA)*HS**(-1.D0/6.D0)*0.4D0*0.5D0*DABS(PS+PP1(J,I))

        AMT(1,JJC) = 0.d0
        
        AMT(2,JJC) = - GRV*DTDY

        AMT(3,JJC) = + DXDY &
                   & + RESIST*DTDXDY

        AMT(4,JJC) = + GRV*DTDY

        AMT(5,JJC) = 0.d0
        
        AMT(6,JJC) = + DXDY*PP1(J,I) &
                   & + GRV*DTDY*(DEP0(J+1,I)-DEP0(J,I))      &
                   & + (PHIN*(PN-PP1(J,I))+PHIS*(PS-PP1(J,I)))*DTDXDDY 
    END DO LPJ1

    JC = MAX0(JC,JJC)
!IF(CD2==1.and.i==i2) THEN
!OPEN(5,FILE = 'AMT.TXT')
!DO J = 1 , JC
!WRITE(5,'(20E20.10)') (AMT(K,J),K=1,6)
!END DO
!CLOSE(5)
!pause
!END IF
    CALL SOLVE5(JC,AMT)
!WRITE(*,*) JC
!IF(CD2==1.and.i==i2) THEN
!OPEN(5,FILE = 'AMT2.TXT')
!DO J = 1 , JC , 1
!WRITE(5,'(20E20.10)') (AMT(K,J),K=1,6)
!END DO
!CLOSE(5)
!PAUSE 'X'
!    END IF

    J0 = 0
    DO J = 2 , JC , 2
        J0 = J0 + 1
        PP2(J0,I) = AMT(6,J)
    END DO
    J0 = 0
    DO J = 1 , JC , 2
        J0 = J0 + 1
        HH2(J0,I) = AMT(6,J)
    END DO
END DO LPI0

RETURN
END SUBROUTINE CONTINUITY_MOMENTUM_Y

