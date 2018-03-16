SUBROUTINE CONTINUITY_MOMENTUM_X(PP1,PP2,HH1,HH2,QQ1,QQ2,T0,DT,DX,DY,I1,I2,IP,UPS)
USE TEST
USE FLOATPRECISION
USE PARAMETERS
USE DEPTH
IMPLICIT NONE
REAL(FP) , DIMENSION(:,:) , POINTER   :: PP1 , PP2 , HH1 , HH2 , QQ1 , QQ2
REAL(FP)                              :: T0  , DT  , DX  , DY
INTEGER                               :: I1  , I2  , IP
LOGICAL                               :: UPS , DWS

REAL(FP) :: FROUDE , WF
REAL(FP) :: ADVE , ADVW
REAL(FP) :: VN , VS , PN , PS , CADL , CADR , CDIFFL , CDIFFR
REAL(FP) :: RESIST , CHEZY , QAV , HAV , GRV , HAV0
REAL(FP) :: HE , HW , HN , HS
REAL(FP) :: DIFFX , DIFFY
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
POW1    = 7.D0/3.D0
COE1    = 4.D0/3.D0

ML = LBOUND(PP1,2) 
MU = UBOUND(PP1,2)
NL = LBOUND(PP1,1) 
NU = UBOUND(PP1,1)


LPI0:DO I = I1 , I2 , IP
    ! CONTINUITY EQUATION
    JC = -1
    LPJ0:DO J = ML , MU , 1
        JC = JC + 2
        IF(NDD(I,J)>0) THEN
            AMT(1,JC) = 0.D0
            AMT(2,JC) = 0.D0
            AMT(3,JC) = 1.D0
            AMT(4,JC) = 0.D0
            AMT(5,JC) = 0.D0
            AMT(6,JC) = HH1(I,J) 
            CYCLE LPJ0   
        END IF
        IF(NDD(I,J)==-2) THEN
            AMT(1,JC) = 0.D0
            AMT(2,JC) = 0.D0
            AMT(3,JC) = 1.D0
            AMT(4,JC) = 0.D0
            AMT(5,JC) = 0.D0
            AMT(6,JC) = HH2(I,J)
            CYCLE LPJ0
        END IF
        HE = 0.5D0*(HH1(I,J)+HH1(I,J+1))
        HW = 0.5D0*(HH1(I,J)+HH1(I,J-1))
        HN = 0.5D0*(HH1(I,J)+HH1(I+1,J))
        HS = 0.5D0*(HH1(I,J)+HH1(I-1,J))
        AMT(1,JC) = 0.D0
        AMT(2,JC) = - DTDY*HW
        AMT(3,JC) = + DXDY4                                   
        AMT(4,JC) = + DTDY*HE
        AMT(5,JC) = 0.D0
        AMT(6,JC) = + DXDY4*HH1(I,J) &                       
                  & - DTDX*(QQ2(I,J)*HN-QQ2(I-1,J)*HS+QQ1(I,J)*HN-QQ1(I-1,J)*HS) &
                  & - DTDY*(PP1(I,J)*HE-PP1(I,J-1)*HW)
    END DO LPJ0

    ! MOMENTUM EQUATION
    JJC = 0
    LPJ1:DO J = ML , MU-1 , 1
        JJC = JJC + 2 ! ID OF BAND MATRIX

        ! CLOSE BOUNDARY CONDITION
        IF(NDD(I,J)>0.OR.NDD(I,J+1)>0.OR.NDDP(I,J)>0) THEN
            AMT(1,JJC) = 0.D0
            AMT(2,JJC) = 0.D0
            AMT(3,JJC) = 1.D0
            AMT(4,JJC) = 0.D0
            AMT(5,JJC) = 0.D0
            AMT(6,JJC) = 0.D0
            CYCLE LPJ1
        END IF

        !  OPEN BOUNDARY CONDITION --- SURFACE LEVEL
        IF(NDD(I,J)==-2.OR.NDD(I,J+1)==-2) THEN
                GRV  = GRA05*(HH1(I,J)+HH1(I,J+1))
                AMT(1,JJC) = 0.D0
                AMT(2,JJC) = - GRV*DTDY
                AMT(3,JJC) = + DXDY*0.5D0*(HH1(I,J)+HH1(I,J+1))
                AMT(4,JJC) = + GRV*DTDY
                AMT(5,JJC) = 0.D0
                AMT(6,JJC) = + DXDY*PP1(I,J)*0.5D0*(HH1(I,J)+HH1(I,J+1)) &
                           & + GRV*DTDY*(DEP0(I,J+1)-DEP0(I,J))  
            CYCLE LPJ1                   
        END IF

        ! FROUDE NUMBER
        FROUDE = DABS(PP1(I,J))/SQRT(GRA0125*(HH1(I,J+1)+HH1(I,J))**3)
        IF(FROUDE>=1.D0) THEN
            WF = 1.D0
        ELSE IF(FROUDE<1.D0.AND.FROUDE>0.25D0) THEN
            WF = COE1*(FROUDE-0.25D0)
        ELSE
            WF = 0.D0
        END IF
        IF(PP1(I,J)<0.D0) WF = -WF

        ! ADVECTION MOMENTUM
        ADVE = 0.25D0*((1.D0-WF)*PP1(I,J+1)+(1.D0+WF)*PP1(I,J)) &
                & / ((1.D0-DABS(WF))*HH1(I,J+1)+DMAX1(WF,0.D0)*HH1(I,J)-DMIN1(WF,0.D0)*HH1(I,J+2))
        ADVW = 0.25D0*((1.D0-WF)*PP1(I,J)+(1.D0+WF)*PP1(I,J-1)) &
                & / ((1.D0-DABS(WF))*HH1(I,J)+DMAX1(WF,0.D0)*HH1(I,J-1)-DMIN1(WF,0.D0)*HH1(I,J+1))
        
        ! CROSS MOMENTUM AND CROSS DIFFUSION
        IF(NDD(I+1,J)==0.AND.NDD(I+1,J+1)==0) THEN
            VN = 2.D0*(QQ2(I,J)+QQ2(I,J+1))/(HH1(I,J)+HH1(I,J+1)+HH1(I+1,J)+HH1(I+1,J+1))
            PN = PP1(I+1,J)
        ELSE IF(NDD(I+1,J)/=0.AND.NDD(I+1,J+1)==0) THEN
            VN = 2.D0*(QQ2(I,J)+QQ2(I,J+1))/(HH1(I,J)+HH1(I,J+1)+HH1(I,J)+HH1(I+1,J+1))
            PN = PP1(I,J)
        ELSE IF(NDD(I+1,J)==0.AND.NDD(I+1,J+1)/=0) THEN
            VN = 2.D0*(QQ2(I,J)+QQ2(I,J+1))/(HH1(I,J)+HH1(I,J+1)+HH1(I+1,J)+HH1(I,J+1))
            PN = PP1(I,J)
        ELSE
            VN = 0.D0
            PN = PP1(I,J)
        END IF
        IF(NDD(I-1,J)==0.AND.NDD(I-1,J+1)==0) THEN
            VS = 2.D0*(QQ2(I-1,J)+QQ2(I-1,J+1))/(HH1(I,J)+HH1(I,J+1)+HH1(I-1,J)+HH1(I-1,J+1))
            PS = PP1(I-1,J)
        ELSE IF(NDD(I-1,J)/=0.AND.NDD(I-1,J+1)==0) THEN
            VS = 2.D0*(QQ2(I-1,J)+QQ2(I-1,J+1))/(HH1(I,J)+HH1(I,J+1)+HH1(I,J)+HH1(I-1,J+1))
            PS = PP1(I,J)
        ELSE IF(NDD(I-1,J)==0.AND.NDD(I-1,J+1)/=0) THEN
            VS = 2.D0*(QQ2(I-1,J)+QQ2(I-1,J+1))/(HH1(I,J)+HH1(I,J+1)+HH1(I-1,J)+HH1(I,J+1))
            PS = PP1(I,J)
        ELSE
            VS = 0.D0
            PS = PP1(I,J)
        END IF
        DIFFY = (VN+VS)**2*DT0125 + EDDYVITY

        IF(.NOT.UPS) THEN
            IF(NDD(I+1,J)<=0.AND.NDD(I+1,J+1)<=0) THEN
                CADL   = - 0.5D0*VS
                CADR   =   0.5D0*(-VN*PP2(I+1,J)+VS*PS-VN*PP1(I,J))
                CDIFFL =   DIFFY
                CDIFFR =   DIFFY*(PP2(I+1,J)-PP1(I,J)+PS)
            ELSE
                CADL   =   0.5D0*(VN-VS)
                CADR   =   0.5D0*(VS*PS-VN*PP1(I,J))        
                CDIFFL =   0.D0
                CDIFFR =   DIFFY*(-PP1(I,J)+PS)
            END IF
        ELSEIF(UPS) THEN
            IF(NDD(I-1,J)<=0.AND.NDD(I-1,J+1)<=0) THEN
                CADL   =   0.5D0*VN
                CADR   =   0.5D0*(-VN*PN+VS*PP2(I-1,J)+VS*PP1(I,J))
                CDIFFL =   DIFFY
                CDIFFR =   DIFFY*(PN-PP1(I,J)+PP2(I-1,J))
            ELSE
                CADL =   0.5D0*(VN-VS)
                CADR =   0.5D0*(-VN*PN+VS*PP1(I,J))
                CDIFFL =   0.D0
                CDIFFR =   DIFFY*(PN-PP1(I,J))     
            END IF
        END IF

        ! GRAVITY 
        GRV  = GRA05*(HH1(I,J+1)+HH1(I,J))
        
        ! RESISTANCE TERM
        QAV = 0.125D0*(QQ1(I-1,J)+QQ1(I-1,J+1)+QQ2(I-1,J)+QQ2(I-1,J+1)+QQ1(I,J)+QQ1(I,J+1)+QQ2(I,J)+QQ2(I,J+1))
        IF(PP1(I,J)>=0.D0) THEN
            HAV = HH1(I,J)
        ELSE IF(PP1(I,J)<0.D0) THEN
            HAV = HH1(I,J+1)
        END IF
        RESIST = GDM2*SQRT(PP1(I,J)**2+QAV**2)/HAV**POW1 

        ! X-DIFFUSION TERM
        DIFFX = (PP1(I,J)/(HH1(I,J+1)+HH1(I,J)))**2*DT2 + EDDYVITY

        AMT(1,JJC) = - (1.D0+WF)*ADVW*DTDY &
                   & - DIFFX*DTDYDDX
        
        AMT(2,JJC) = - GRV*DTDY

        AMT(3,JJC) = + DXDY &
                   & + ((1.D0+WF)*ADVE-(1.D0-WF)*ADVW)*DTDY &
                   & + CADL*DTDX &
                   & + RESIST*DTDXDY &
                   & + 2.D0*DIFFX*DTDYDDX &
                   & + CDIFFL*DTDXDDY

        AMT(4,JJC) = + GRV*DTDY

        AMT(5,JJC) = + (1.D0-WF)*ADVE*DTDY &
                   & - DIFFX*DTDYDDX
        
        AMT(6,JJC) = + DXDY*PP1(I,J) &
                   & + CADR*DTDX &
                   & + CDIFFR*DTDXDDY  &
                   & + GRV*DTDY*(DEP0(I,J+1)-DEP0(I,J))
        
    END DO LPJ1

    JC = MAX0(JC,JJC)
!    IF(I==59) THEN
!OPEN(5,FILE = 'AMT.TXT')
!DO J = 1 , JC
!WRITE(5,'(20E20.10)') (AMT(K,J),K=1,6)
!END DO
!CLOSE(5)
!WRITE(*,*) I
!    END IF
    CALL SOLVE5(JC,AMT)
    
!IF(I==59) THEN
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
        PP2(I,J0) = AMT(6,J)
    END DO
    J0 = 0
    DO J = 1 , JC , 2
        J0 = J0 + 1
        HH2(I,J0) = AMT(6,J)
    END DO
END DO LPI0

RETURN
END SUBROUTINE CONTINUITY_MOMENTUM_X

