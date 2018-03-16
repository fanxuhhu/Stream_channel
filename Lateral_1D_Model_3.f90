#include "AMYCAI.h"

MODULE Lateral_1D_Model

    USE TimeCatcher
    USE FloatPrecision
    USE ArraySize
    USE GRIDS
    USE hydro
    USE ReGrid
    USE Parameters
    USE LinearSolver
    
    IMPLICIT NONE
    
    CONTAINS

!= = = = = = = = = = = = = = = = = = = = = = = = = = 
!
!   Main Routine to calculate the transversal 
!   distribution of shear stress and velocity
!
!= = = = = = = = = = = = = = = = = = = = = = = = = = 

    SUBROUTINE Cross_Sectional_Tau(TDUD,TDTAU,QQ1,TDEP0,Upsilon,XDistance,WetSegments,FC,MaxTau,MeanTau,Sf,RHOGSD)

        IMPLICIT NONE
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: TDUD
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: TDTAU
        REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(IN)    :: QQ1
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(IN)    :: TDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) , TARGET  , INTENT(IN)    :: Upsilon
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,           INTENT(IN)    :: XDistance
        INTEGER  , ALLOCATABLE , DIMENSION(:)   ,           INTENT(INOUT) :: WetSegments
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET  , INTENT(INOUT) :: FC      
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)   ,         INTENT(INOUT) :: MaxTau
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)   ,         INTENT(INOUT) :: MeanTau
        REAL(fp) , ALLOCATABLE , DIMENSION(:)     ,         INTENT(INOUT) :: Sf
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)   ,         INTENT(INOUT) :: RHOGSD
        !--------------------------------------------------------------------------------  
        
        INTEGER                                                           :: Re_Multiples 
        INTEGER                                                           :: NG
        INTEGER                                                           :: I   , J  , JJ , J1 , J2
        REAL(fp)                                                          :: ReX , XL , XR , Dep1
        INTEGER                                                           :: status
        CHARACTER(len=80)                                                 :: err_msg

!- - - - - - - - - - - - - - - - - - - - - - - - - -  
!
!   Calculate the shear stress over each cross-section
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - 

        DO I = ML , MU

! - - - - - - - - - - - - - - - - - - - -     
! get the wet points of each cross-section
! - - - - - - - - - - - - - - - - - - - -

            !JSTART(I) = NL
            !JEND(I) = NU
            !LPJ9: DO J = (NL+NU)/2 , NU
            !    IF(TDEP0(I,J)<0.D0) THEN
            !        JEND(I)  = J
            !        JSTART(I) = NL + NU - JEND(I)
            !        EXIT LPJ9
            !    END IF          
            !END DO LPJ9
            !IF(JSTART(I)>(NU+NL)/2) THEN
            !    PAUSE 'error regrid'
            !END IF    
            !WetSegments(I) = -JSTART(I) + JEND(I) - 1

            JSTART(I) = NL
            JEND(I) = NU
            LPJ9: DO J = NU , NL , -1
                IF(TDEP0(I,J)<0.D0) THEN
                    JSTART(I) = J
                    EXIT LPJ9
                END IF          
            END DO LPJ9
            WetSegments(I) = -JSTART(I) + JEND(I) - 1
            
! - - - - - - - - - - - - - - - - - - - -     
! IF the WetSegments > 9, directly calculate
! the shear stress
! - - - - - - - - - - - - - - - - - - - -   
            
            IF(WetSegments(I)>=9) THEN

                IF(QQ1(I)/=0.D0) THEN
                    CALL Iterate_Tau(I,NL,NU,DX,QQ1(I),TDEP0(I,:),Upsilon(I,:),TDUD(I,:),TDTAU(I,:),FC(I),Sf(I))
                ELSEIF(QQ1(I)==0.D0) THEN
                    TDUD(I,:)     = 0.D0
                    TDTAU(I,:) = 0.D0
                    FC(I)     = 0.D0
                    Sf(I) = 0.D0 
                END IF
                
! - - - - - - - - - - - - - - - - - - - -     
! IF the WetSegments < 9, refine the grid
! - - - - - - - - - - - - - - - - - - - -   

            ELSEIF(WetSegments(I)<9.and.WetSegments(I)>0) THEN

! calculate multiple times
                
                Re_Multiples = 9 / WetSegments(I) + 1
                
! number of new points

                NG = Re_Multiples * (WetSegments(I) + 1) + 1

! New Grid size
            
                Re_dx = dx / DBLE(Re_Multiples)
                
! Calculate TTDepthQ after regrid

                !IF(WetSegments(I)>5) THEN
            
                    DO J = 1 , NG
                        ReX = XDistance(JSTART(I)) + dble(J-1)*Re_dx      
                        LPR1:DO JJ = JSTART(I) , JEND(I)-1
                            XL = XDistance(JJ)
                            XR = XDistance(JJ+1)
                            IF(ReX>=XL.and.ReX<=XR) THEN
                                Re_TTDep(J) = ((ReX-XL)*TDEP0(I,JJ+1)+(XR-ReX)*TDEP0(I,JJ))/dx
                                Re_Upsilon0(J) = ((ReX-XL)*Upsilon(I,JJ+1)+(XR-ReX)*Upsilon(I,JJ))/dx
                                EXIT LPR1
                            END IF
                        END DO LPR1
                    END DO

                !ELSE
                !
                !    Re_Upsilon0(:) = 1.D0
                !    Re_TTDep(:) = TDEP0(I,NU)
                !    Re_TTDep(1) = TDEP0(I,JSTART(I))
                !    Re_TTDep(NG) = TDEP0(I,JEND(I))
                !
                !ENDIF
            
                IF(QQ1(I)/=0.D0) THEN
                
                    CALL Iterate_Tau(I,NL,NG,Re_dx,QQ1(I),Re_TTDep,Re_Upsilon0,Re_U0,Re_TAU0,FC(I),Sf(I))
                
! fit back to coarse grid
                
                    JJ = JSTART(I) - 1
                    DO J = 1 , NG , Re_Multiples
                        JJ = JJ + 1
                        TDUD(I,JJ) = Re_U0(J)
                        TDTAU(I,JJ) = Re_TAU0(J)
                    END DO
                
                ELSEIF(QQ1(I)==0.D0) THEN
                    
                    TDUD(I,:)   = 0.D0
                    TDTAU(I,:) = 0.D0
                    FC(I)     = 0.D0
                    Sf(I) = 0.D0
                    
                END IF

            ELSE
                
                TDUD(I,:)   = 0.D0
                TDTAU(I,:) = 0.D0
                FC(I)     = 0.D0                
                Sf(I) = 0.D0    
                    
            ENDIF
                    
        END DO
 
! MaxTau

        DO I = ML , MU
            DO J = NL , NU
                MaxTau(I,J)  = dmax1(MaxTau(I,J),TDTAU(I,J))
                MeanTau(I,J) = MeanTau(I,J) + TDTAU(I,J)
            ENDDO
        ENDDO
       
!RHOGSD
        DO I = ML , MU
            DO J = NL , NU
                IF(J==NU) THEN
                    J1 = J
                ELSE
                    J1 = J+1
                ENDIF
                
                IF(J==NL) THEN
                    J2 = J
                ELSE
                    J2 = J-1
                ENDIF            
            
                IF(TDEP0(I,J)>0.D0) THEN
                    Dep1 = TDEP0(I,J)
                ELSE
                    Dep1 = DMAX1(TDEP0(I,J1),TDEP0(I,J2))
                ENDIF   
                
                RHOGSD(I,J) = RHOW*gra*Dep1*DABS(SF(I))
                
            ENDDO
        ENDDO
        
        
END SUBROUTINE Cross_Sectional_Tau


!= = = = = = = = = = = = = = = = = = = = = = = = = = 
!
!   Iterate the shear stress over each cross-section
!
!= = = = = = = = = = = = = = = = = = = = = = = = = = 


    SUBROUTINE Iterate_Tau(IROW,JLOW,JUP,LX0,TargetQ,TTDep,Upsilon0,U0,TAU0,FC0,SFS)

        IMPLICIT NONE

        INTEGER                                         , INTENT(IN)    :: IROW   
        INTEGER                                         , INTENT(IN)    :: JLOW
        INTEGER                                         , INTENT(IN)    :: JUP
        REAL(fp)                                        , INTENT(IN)    :: LX0
        REAL(fp)                                        , INTENT(IN)    :: TargetQ        
        REAL(fp)               , DIMENSION(:)           , INTENT(IN)    :: TTDep
        REAL(fp)               , DIMENSION(:)           , INTENT(IN)    :: Upsilon0
        REAL(fp)               , DIMENSION(:)           , INTENT(INOUT) :: U0
        REAL(fp)               , DIMENSION(:)           , INTENT(INOUT) :: TAU0
        REAL(fp)                                        , INTENT(INOUT) :: FC0   
        REAL(fp)                                        , INTENT(INOUT) :: SFS
        !-----------------------------------------------------------------------------
        REAL(fp)               , PARAMETER                              :: TOL = 1.D-6
        REAL(fp)                                                        :: ShootQ       
        REAL(fp)                                                        :: Dep1 , FF1  
        REAL(fp)                                                        :: CC1  , CC2
        INTEGER                                                         :: J,ITE,JC,JJ,J1,J2     
        INTEGER                                                         :: status
        CHARACTER(LEN=80)                                               :: err_msg
        REAL(fp) ,  DIMENSION(500) :: LateralDiff
        REAL(fp) ,  DIMENSION(500) :: DWFriction

! IF discharge is zero, all is zero

        IF(TargetQ==0) THEN
            RETURN
        END IF   
        
! Calculate the friction coefficients
        
        DO J = JLOW , JUP
            
            IF(J==JUP) THEN
                J1 = J
            ELSE
                J1 = J+1
            ENDIF
                
            IF(J==JLOW) THEN
                J2 = J
            ELSE
                J2 = J-1
            ENDIF            
            
            IF(TTDep(J)>0.D0) THEN
                Dep1 = TTDep(J)
            ELSE
                Dep1 = DMAX1(TTDep(J1),TTDep(J2))
            ENDIF            
            
            IF(Dep1>0.D0) THEN
#ifdef Shiono_Knight_1988
                !DWFriction(J) = manning**2*gra*DMAX1(TTDep(J),0.1D0)**(-1.D0/3.D0)
                DWFriction(J) = manning**2*gra*Dep1**(-1.D0/3.D0) !* Upsilon0(J)
#else
                DWFriction(J) = manning**2*gra*Dep1**(-1.D0/3.D0)
#endif
            ELSE
                DWFriction(J) = 0.D0
            END IF
        END DO

! Calculate the diffusion term
        
        DO J = JLOW , JUP-1
            IF(TTDep(J)>0.D0.AND.TTDep(J+1)>0.D0) THEN
                !Dep1 = DMAX1(0.5d0*(TTDep(J)+TTDep(J+1)),0.5D0)
                Dep1 = 0.5d0*(TTDep(J)+TTDep(J+1))
                FF1 = manning**2*gra*Dep1**(-1.D0/3.D0)
                LateralDiff(J) = 0.5D0*Lambda*Dep1**2*SQRT(FF1)
            ELSEIF(TTDep(J)>0.D0.AND.TTDep(J+1)<=0.D0) THEN
                !Dep1 = DMAX1(TTDep(J),0.5D0)
                Dep1 = TTDep(J)
                FF1 = manning**2*gra*Dep1**(-1.D0/3.D0)
                LateralDiff(J) = 0.5D0*Lambda*Dep1**2*SQRT(FF1)
            ELSEIF(TTDep(J+1)>0.D0.AND.TTDep(J)<=0.D0) THEN
                !Dep1 = DMAX1(TTDep(J+1),0.5D0)
                Dep1 = TTDep(J+1)
                FF1 = manning**2*gra*Dep1**(-1.D0/3.D0)
                LateralDiff(J) = 0.5D0*Lambda*Dep1**2*SQRT(FF1)
            ELSE
                LateralDiff(J) = 0.D0
            END IF
        END DO

! initial guess of SFS : Surface slope
        
        SFS = 0.5D0

!  - - - - - - - - - - - - - - - - - - - - - -
! LPK: the iteration loop
! The upper limit of iteration times is 1000
!- - - - - - - - - - - - - - - - - - - - - - -
        
        LPK: DO ITE = 1 , 1000
 
            JC = 0

! - - - - - - - - - - - - - - - - - - - - - - - - 
! LPJ: the loop calculate the transversal distribution 
!      of the shear stress AND the velocity
! - - - - - - - - - - - - - - - - - - - - - - - -        
    
            LPJ:DO J = JLOW , JUP
                
                JC = JC + 1

                IF(J==JUP) THEN
                    J1 = J
                ELSE
                    J1 = J+1
                ENDIF
                
                IF(J==JLOW) THEN
                    J2 = J
                ELSE
                    J2 = J-1
                ENDIF
                
! Dry point
                
                IF(TTDep(J)<=0.AND.TTDep(J1)<=0.AND.TTDep(J2)<=0) THEN
                    AMT(1,JC) = 0.D0
                    AMT(2,JC) = 1.D0
                    AMT(3,JC) = 0.D0
                    AMT(4,JC) = 0.D0
                    cycle LPJ
                END IF

! DEPTH
                
                IF(TTDep(J)>0.D0) THEN
                    Dep1 = TTDep(J)
                ELSE
                    Dep1 = DMAX1(TTDep(J1),TTDep(J2))
                ENDIF
                
! Left Edge (Intertidal side)    
                
                IF(J==JLOW) THEN
                    CC2 = 1.D0/LX0**2*LateralDiff(J)
                    AMT(1,JC) = 0.D0       
                    AMT(2,JC) = DWFriction(J) + 2.D0*CC2
                    AMT(3,JC) = - 2.D0*CC2
                    AMT(4,JC) = gra*Dep1*DABS(SFS)
                    cycle LPJ
                END IF
                
! Right Edge (Intertidal side)
                
                IF(J==JUP) THEN
                    CC1 = 1.D0/LX0**2*LateralDiff(J-1)
                    AMT(1,JC) = - 2.D0*CC1       
                    AMT(2,JC) = DWFriction(J) + 2.D0*CC1
                    AMT(3,JC) = 0.D0
                    AMT(4,JC) = gra*Dep1*DABS(SFS)
                    cycle LPJ
                END IF
                
! General case
                
                CC1 = 1.D0/LX0**2*LateralDiff(J-1) ! DIFFUSION TERMS  
                CC2 = 1.D0/LX0**2*LateralDiff(J)

                IF(CC1==0.D0.AND.CC2>0.D0) THEN ! Left dry
                
                    AMT(1,JC) = 0.D0   
                    AMT(2,JC) = DWFriction(J) + 2.D0*CC2 
                    AMT(3,JC) = - 2.D0*CC2 
                    AMT(4,JC) = gra*Dep1*DABS(SFS)
                    
                ELSEIF(CC1>0.D0.AND.CC2==0.D0) THEN ! Right dry
                    
                    AMT(1,JC) = - 2.D0*CC1
                    AMT(2,JC) = DWFriction(J) + 2.D0*CC1 
                    AMT(3,JC) = 0.D0 
                    AMT(4,JC) = gra*Dep1*DABS(SFS)  
                    
   
                ELSE ! All wet 
                    
                    AMT(1,JC) = - CC1        
                    AMT(2,JC) = DWFriction(J) + CC1 + CC2 
                    AMT(3,JC) = - CC2
                    AMT(4,JC) = gra*Dep1*DABS(SFS)
                
                END IF                

            END DO LPJ

            !IF(IROW==1.and.cd==1021162) then
            !OPEN(111,file='AMT1.txt')
            !DO J = 1 , JC
            !    WRITE(111,'(50e25.15)') (AMT(JJ,J),JJ=1,4),DWFriction(J)
            !END DO
            !CLOSE(111)
            !END IF

            CALL solve3_LU(JC)
            
            !IF(IROW==1.and.cd==1021162) then
            !OPEN(111,file='AMT2.txt')
            !DO J = 1 , JC
            !    WRITE(111,'(50e25.15)') (AMT(JJ,J),JJ=1,4)
            !END DO
            !CLOSE(111)
            !pause 'AMT6'
            !END IF
            
! get the velocity

            DO J = JLOW , JUP
                U0(J) = SQRT(AMT(4,J))
                IF(TargetQ<0.D0) U0(J) = - U0(J)
            END DO

! estimate the discharge

            ShootQ = 0.D0
            DO J = JLOW , JUP-1
                ShootQ = ShootQ + 0.25D0*(U0(J)+U0(J+1))*(dmax1(0.D0,TTDep(J+1))+dmax1(0.D0,TTDep(J)))*LX0
            END DO

! compare with the actual discharge

            IF(DABS(ShootQ-TargetQ) <= TOL) THEN

                DO J = 1 , JC
                    TAU0(J) = DWFriction(J)*U0(J)**2*RHOW !/ Upsilon0(J)                
                END DO
                    
! calculate FC
! Coupling the lateral 1d model with the longitudinal model.
! I.e., modify the friction coefficient of the longitudinal
! 1d model to match the intergral bed resistance calculated
! by 1d lateral model.
                    
                !FC0 = 0.D0
                !
                !DO J = JLOW , JUP-1 , 1
                !    IF(DWFriction(J)>0.D0.OR.DWFriction(J+1)>0.D0) THEN
                !        FC0 = FC0 &
                !                & + 0.5D0 * (&
                !                &   DWFriction(J)*U0(J)**2 &
                !                & + DWFriction(J+1)*U0(J+1)**2 &
                !                & ) * LX0
                !    END IF
                !END DO
                !FC0 = FC0 / DABS(TargetQ)         
                
                RETURN

            END IF

            SFS = SFS*(TargetQ/ShootQ)**2 ! modify the energy slope

        END DO LPK

! error message

        IF(ITE == 1001) THEN
            WRITE(*,*) IROW , CD
            PAUSE 'error1'
        END IF

    END SUBROUTINE Iterate_Tau

END MODULE Lateral_1D_Model