SUBROUTINE CONTINUITY_MOMENTUM_1D(QQ1,QQ2,EE1,EE2,AA1,BB1,GAMA,T0,DT,DX,ML,MU)
    use TimeCatcher
    use FloatPrecision
    use Parameters
    use Topography
    implicit none

    integer                     :: ML , MU    
    real(fp) , dimension(ML:MU) :: QQ1 , QQ2 , EE1 , EE2 , AA1 , BB1 , GAMA
    real(fp)                    :: T0  , DT  , DX  , DY

    
    real(fp) :: FROUDE , WF
    real(fp) :: ADVE , ADVW
    real(fp) :: RESIST, AVA , GRV
    real(fp) :: DIFFX
    integer  :: I , J , K , JC , JJC , J0
    real(fp) :: AMT(6,3001)
    real(fp) :: BndSeaLevel

! = = = = = = = = = = = = = = = = = = = = = = =
!
!           CONTINUITY EQUATION LOOP: LPJ0
!
! = = = = = = = = = = = = = = = = = = = = = = =
    
    JC = -1
    LPJ0:DO J = ML , MU , 1
        JC = JC + 2
        IF(ND(J)==-2) THEN
            AMT(1,JC) = 0.D0
            AMT(2,JC) = 0.D0
            AMT(3,JC) = 1.D0
            AMT(4,JC) = 0.D0
            AMT(5,JC) = 0.D0
            AMT(6,JC) = BndSeaLevel(T0,TDAMP,TDFRQ)
            CYCLE LPJ0
        end IF
        IF(ND(J)>0.or.(NDQ(J)>0.AND.NDQ(J-1)>0).or.(J==MU.AND.NDQ(J-1)>0)) THEN
            AMT(1,JC) = 0.D0
            AMT(2,JC) = 0.D0
            AMT(3,JC) = 1.D0
            AMT(4,JC) = 0.D0
            AMT(5,JC) = 0.D0
            AMT(6,JC) = EE1(J) 
            CYCLE LPJ0   
        end IF        
        AMT(1,JC) =   0.D0
        AMT(2,JC) = - DT
        AMT(3,JC) = + DX*BB1(J)                                   
        AMT(4,JC) = + DT
        AMT(5,JC) =   0.D0
        AMT(6,JC) = + DX*BB1(J)*EE1(J) &                       
                  & - DT*(QQ1(J)-QQ1(J-1))
    end DO LPJ0

! = = = = = = = = = = = = = = = = = = = = = = =
!    
!        MOMENTUM EQUATION LOOP: LPJ1
!    
! = = = = = = = = = = = = = = = = = = = = = = =
    
    JJC = 0
    LPJ1:DO J = ML , MU-1 , 1
        
        JJC = JJC + 2
        
!----------------------------------------------
!        Close boundary condition
!----------------------------------------------
        
        IF(ND(J)>0.OR.ND(J+1)>0.OR.NDQ(J)>0) THEN
            AMT(1,JJC) = 0.D0
            AMT(2,JJC) = 0.D0
            AMT(3,JJC) = 1.D0
            AMT(4,JJC) = 0.D0
            AMT(5,JJC) = 0.D0
            AMT(6,JJC) = 0.D0
            CYCLE LPJ1
        end IF

!----------------------------------------------
!   Open boundary condition -- surface level
!----------------------------------------------
        
        IF(ND(J)==-2.OR.ND(J+1)==-2) THEN
            GRV = GRA*0.5D0*(AA1(J)+AA1(J+1))
            AMT(1,JJC) = 0.D0
            AMT(2,JJC) = - GRV*DT
            AMT(3,JJC) = + DX + GAMA(J)
            AMT(4,JJC) = + GRV*DT
            AMT(5,JJC) = 0.D0
            AMT(6,JJC) = + DX*QQ1(J)
            CYCLE LPJ1                   
        end IF

!----------------------------------------------
!                General Case
!----------------------------------------------

! Averaged cross-sectional area
    
        AVA = 0.5D0*(AA1(J+1)+AA1(J))

! Froude Number
        
        FROUDE = DABS(QQ1(J))*SQRT(BB1(J)/GRA/AVA)
        IF(FROUDE>=1.D0) THEN
            WF = 1.D0
        ELSE IF(FROUDE<1.D0.AND.FROUDE>0.25D0) THEN
            WF = (3.D0/4.D0)*(FROUDE-0.25D0)
        ELSE
            WF = 0.D0
        end IF
        IF(QQ1(J)<0.D0) WF = -WF

! Advection term

    ADVE = 0.25D0*((1.D0-WF)*QQ1(J+1)+(1.D0+WF)*QQ1(J)) &
            & / ((1.D0-DABS(WF))*AA1(J+1)+DMAX1(WF,0.D0)*AA1(J)-DMIN1(WF,0.D0)*AA1(J+2))
    ADVW = 0.25D0*((1.D0-WF)*QQ1(J)+(1.D0+WF)*QQ1(J-1)) &
            & / ((1.D0-DABS(WF))*AA1(J)+DMAX1(WF,0.D0)*AA1(J-1)-DMIN1(WF,0.D0)*AA1(J+1))
        
! Gravity term
        
        GRV  = GRA*AVA
        
! Resistance term
        
        RESIST = GAMA(J)

! X-Diffusion term

        DIFFX = 0.5D0*(QQ1(J)/AVA)**2*DT
        
! Fill the matrix    
        
        AMT(1,JJC) = - (1.D0+WF)*ADVW*DT &
                   & - DIFFX*DT/DX
        
        AMT(2,JJC) = - GRV*DT

        AMT(3,JJC) = + DX &
                   & + ((1.D0+WF)*ADVE-(1.D0-WF)*ADVW)*DT &
                   & + RESIST*DT*DX &
                   & + 2.D0*DIFFX*DT/DX

        AMT(4,JJC) = + GRV*DT

        AMT(5,JJC) = + (1.D0-WF)*ADVE*DT &
                   & - DIFFX*DT/DX
        
        AMT(6,JJC) = + DX*QQ1(J)
        
    end DO LPJ1

    JC = MAX0(JC,JJC)
!    IF(I==59) THEN
!OPEN(5,FILE = 'AMT3.TXT')
!DO J = 1 , JC
!WRITE(5,'(20E20.10)') (AMT(K,J),K=1,6)
!end DO
!CLOSE(5)
!WRITE(*,*) I
!    end IF
    CALL SOLVE5(JC,AMT)
    
!IF(I==59) THEN
!OPEN(5,FILE = 'AMT4.TXT')
!DO J = 1 , JC , 1
!WRITE(5,'(20E20.10)') (AMT(K,J),K=1,6)
!end DO
!CLOSE(5)
!PAUSE 'X'
!    end IF

    J0 = 0
    DO J = 2 , JC , 2
        J0 = J0 + 1
        QQ2(J0) = AMT(6,J)
    end DO
    J0 = 0
    DO J = 1 , JC , 2
        J0 = J0 + 1
        EE2(J0) = AMT(6,J)
    end DO

end SUBROUTINE CONTINUITY_MOMENTUM_1D

