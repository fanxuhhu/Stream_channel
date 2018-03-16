! This MODULE contains all the routines that 
! are related to sediment transport. Only the diffusion
! proccess is considered for the 1st verson.    
    
#include "AMYCAI.h"
    
MODULE Sediment_Transport

    USE TimeCatcher
    USE FloatPrecision
    USE ArraySize
    USE GRIDS
    USE hydro
    USE Parameters
    USE Suspended_Load
    USE LinearSolver
    
    IMPLICIT NONE
    
    CONTAINS

! = = = = = = = = = = = = = = = =    
!
!      1D diffusion model
!
! = = = = = = = = = = = = = = = =

    SUBROUTINE MUD_Transport(TDUD,TDTAU,MAXTAU,CRTAU,STDEP,HND,HNDD,HDIFF,HMTDEP0,Q_RESIDUAL1,SSC1,SSC2,SSC1D1,SSC1D2,ERO,TTERO,ERO1D,LABL,LAAVA,EVOSTEP,DEP0,UPSILON,BC,Sf,T0,TDEP0)
    
        IMPLICIT NONE
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: TDUD       
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: TDTAU
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: MAXTAU
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:,:) ,         INTENT(INOUT) :: CRTAU 
        REAL(fp) , ALLOCATABLE , DIMENSION(:) ,             INTENT(IN)    :: STDEP 
        INTEGER  , ALLOCATABLE , DIMENSION(:)   ,           INTENT(IN)    :: HND
        INTEGER  , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: HNDD
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,           INTENT(IN)    :: HDIFF
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,           INTENT(IN)    :: HMTDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,           INTENT(IN)    :: Q_RESIDUAL1 
        REAL(fp) ,               DIMENSION(:,:) , POINTER , INTENT(INOUT) :: SSC1
        REAL(fp) ,               DIMENSION(:,:) , POINTER , INTENT(INOUT) :: SSC2
        REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: SSC1D1
        REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(INOUT) :: SSC1D2
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: ERO
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: TTERO
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,           INTENT(INOUT) :: ERO1D
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: LABL
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: LAAVA 
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: EVOSTEP
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: DEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) , TARGET  , INTENT(IN)    :: UPSILON
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(IN)    :: BC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)             , INTENT(INOUT) :: Sf
        REAL(fp)                                          , INTENT(IN)    :: T0       
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(INOUT) :: TDEP0
        !--------------------------------------------------------------------------------
        REAL(fp) ,               DIMENSION(:)   , POINTER                 :: TEMP

!- - - - - - - - - - - - - -
! Calculate local erosion        
!- - - - - - - - - - - - - -
        CALL LocalErosion(ERO,TTERO,ERO1D,TDUD,TDTAU,UPSILON,CRTAU,STDEP,DEP0)
!- - - - - - - - - - - - - -
! Calculate local erosion        
!- - - - - - - - - - - - - -
#ifdef SSC_CONSTANT
        CALL LocalSSC(Sf,DEP0,SSC1,SSC2,T0,TDEP0)
        !CALL SSC_MUDFLAT(Sf,DEP0,SSC1,SSC2)
#endif
! - - - - - - - - - - - - - 
!   Calculate the Suspended sediment transport
!- - - - - - - - - - - - - -        
#ifdef SSC_DIFFUSION_1D
        !! 1D solver            
        !CALL SSC1D_Solver(ML,MU,DT,DX,SSC1D1,SSC1D2,HND,HDIFF,ERO1D,HMTDEP0,Q_RESIDUAL1)
        !! Exchange pointers       
        !TEMP=>SSC1D2
        !SSC1D2=>SSC1D1
        !SSC1D1=>TEMP        
        !SSC1D1(:) = ERO1D(:) / WS
        !CALL SSC1D_DIFF_0ORDER(ML,MU,DT,DX,SSC1D1,SSC1D2,HND,100.D0)
#endif
!- - - - - - - - - - - - - - - - 
! Lateral adjustment
!- - - - - - - - - - - - - - - -
    ! Bed load Transport
#if defined(BEDLOAD_INC) .OR. defined(AVAlANCHE_INC)
        CALL LateralST(DEP0,TDTAU,LABL,LAAVA)
#endif
!- - - - - - - - - - - - - - - - 
! bed evolution
!- - - - - - - - - - - - - - - - 
#ifdef LONGTERM
        !IF(CD2>=15) THEN
#ifdef MAXTAUEVO
            CALL BedEvolution(MAXTAU,ERO,SSC1D1,LABL,LAAVA,HNDD,EVOSTEP,DEP0)
#else
            CALL BedEvolution_im(TDTAU,CRTAU,STDEP,ERO,SSC1,SSC1D1,LABL,LAAVA,HNDD,EVOSTEP,DEP0)
#endif
       ! ENDIF
#endif

    END SUBROUTINE MUD_Transport

!= = = = = = = = 
!
! 1D sediment diffusion equation
!
!= = = = = = = =  

    SUBROUTINE SSC1D_Solver(IL,IU,LT,LX,CP1,CP2,GP,DF,EP,HP,QP)
    
        IMPLICIT NONE
                                                        
        INTEGER  ,                                        INTENT(IN)    :: IL
        INTEGER  ,                                        INTENT(IN)    :: IU
        REAL(fp) ,                                        INTENT(IN)    :: LT
        REAL(fp) ,                                        INTENT(IN)    :: LX
        REAL(fp) ,               DIMENSION(:) , POINTER , INTENT(IN)    :: CP1
        REAL(fp) ,               DIMENSION(:) , POINTER , INTENT(INOUT) :: CP2
        INTEGER  , ALLOCATABLE , DIMENSION(:) ,           INTENT(IN)    :: GP
        REAL(fp) , ALLOCATABLE , DIMENSION(:) ,           INTENT(IN)    :: DF
        REAL(fp) , ALLOCATABLE , DIMENSION(:) ,           INTENT(IN)    :: EP
        REAL(fp) , ALLOCATABLE , DIMENSION(:) ,           INTENT(IN)    :: HP
        REAL(fp) , ALLOCATABLE , DIMENSION(:) ,           INTENT(IN)    :: QP     
        !----------------------------------------------------------------------
        REAL(fp)                                                        :: LXX        
        REAL(fp)                                                        :: DFE , DFW
        REAL(fp)                                                        :: ADE , ADW , ADP
        INTEGER                                                         :: GDE, GDW
        INTEGER                                                         :: I , IC , K
        
        IC = 0
        LXX = LX**2    
        
        DO I = IL , IU
            
            IC = IC+1
! Local point is dry, set zero
            IF(GP(I)==1) THEN
                amt(1,IC) = 0.D0
                amt(2,IC) = 1.D0
                amt(3,IC) = 0.D0
                amt(4,IC) = 0.D0
                CYCLE      
! Local point is wet, set diffusion coefficients
            ELSEIF(GP(I)==0) THEN
                !WEST
                IF(I==IL) THEN
                    GDW = 1    ! NEUMANN BOUNDARY
                    !GDW = -2    ! zero Dirichlet
                ELSE
                    GDW = GP(I-1)
                END IF
                IF(GDW==0) THEN
                    DFW = DF(I-1) * LT / LXX / HP(I)  
                END IF                  
                !EAST                
                IF(I==IU) THEN
                    GDE = 1
                ELSE
                    GDE = GP(I+1)
                END IF                                
                IF(GDE==0) THEN
                    DFE = DF(I)  * LT / LXX / HP(I)
                END IF
!- - - - - - - - - - - - - - - - 
! Residual Advection --- Upwind scheme
!- - - - - - - - - - - - - - - - 
                IF(I/=IU.and.I/=IL) THEN
                    ADE = DMIN1(QP(I),0.D0) * LT / LX 
                    ADW = -DMAX1(QP(I-1),0.D0) * LT / LX
                    ADP = (DMAX1(QP(I),0.D0)-DMIN1(QP(I-1),0.D0)) * LT / LX
                ELSE
                    ADE = 0.D0
                    ADW = 0.D0
                    ADP = 0.D0
                END IF
!- - - - - - - - - - - - - - - - 
! Left hand side of the equation
!- - - - - - - - - - - - - - - -                 
                IF(GDE==1.AND.GDW==0) THEN                    
                    amt(1,IC) = - 2.D0*DFW
                    amt(2,IC) = + 2.D0*DFW + 1.D0
                    amt(3,IC) = 0.D0                   
                ELSEIF(GDW==1.AND.GDE==0) THEN
                    amt(1,IC) =   0.D0
                    amt(2,IC) = + 2.D0*DFE + 1.D0
                    amt(3,IC) = - 2.D0*DFE                    
                ELSEIF(GDW==1.AND.GDE==1) THEN
                    amt(1,IC) =   0.D0
                    amt(2,IC) =   1.D0
                    amt(3,IC) =   0.D0
                ELSEIF(GDW==0.AND.GDE==0)  THEN
                    amt(1,IC) = - DFW
                    amt(2,IC) = + DFW + DFE + 1.D0
                    amt(3,IC) = - DFE
                ELSEIF(GDW==-2.AND.GDE==0)  THEN   ! zero dirichlet boundary
                    amt(1,IC) =   0.D0
                    amt(2,IC) = + DFW + DFE + 1.D0
                    amt(3,IC) = - DFE    
                ELSE
                    WRITE(*,*) I , GP(I-1) , GP(I+1)
                    PAUSE 'tri_matrix_error_1'
                END IF
!- - - - - - - - - - - - - - - - 
! Advection term ---- RESIDUAL
!- - - - - - - - - - - - - - - - 
#ifdef RESIDUALINC  
                !amt(1,IC) = amt(1,IC) + ADW
                !amt(2,IC) = amt(2,IC) + ADP
                !amt(3,IC) = amt(3,IC) + ADE
#endif
!- - - - - - - - - - - - - - - - 
! Right hand side of the equation
!- - - - - - - - - - - - - - - - 
                amt(4,IC) = CP1(I)
        ! sourse term
                amt(4,IC) = amt(4,IC) + EP(I) / HP(I) * LT
                amt(2,IC) = amt(2,IC) + WS    / HP(I) * LT
            END IF
            
        END DO

        !if(cd==50) then
        !OPEN(5,file='amt3.txt')
        !do I = 1, IC
        !    WRITE(5,'(10e25.10)') (amt(K,I),K=1,4)
        !end do
        !CLOSE(5)
        !end if
        
        CALL solve3_LU(IC)

        !if(cd==50) then
        !OPEN(5,file='amt4.txt')
        !do I = 1, IC
        !    WRITE(5,'(10e25.10)') (amt(K,I),K=1,4)
        !end do
        !CLOSE(5)
        !pause
        !end if
        
        CP2(:) = amt(4,:)

    END SUBROUTINE SSC1D_Solver    

! - - - - - - - - - - - - - - - - 
!
!       Local Erosion
!
! - - - - - - - - - - - - - - - - 

    SUBROUTINE LocalErosion(ERO,TTERO,ERO1D,TDUD,TDTAU,UPSILON,CRTAU,STDEP,DEP0)
    
        IMPLICIT NONE

        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: ERO
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: TTERO
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,           INTENT(INOUT) :: ERO1D        
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: TDUD
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: TDTAU
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) , TARGET  , INTENT(IN)    :: Upsilon
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:,:) ,         INTENT(INOUT) :: CRTAU 
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,           INTENT(IN)    :: STDEP
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: DEP0 
        
        INTEGER :: I , J , K
        REAL(fp) :: taue0 , h0 , erosion0 , erosion1 , tt0 , tt1

    ! tau >  critical shear stress: Linear Function
    ! tau <= critical shear stress: No Erosion
                
#ifdef STRATIFICATION

        tt0 = MorFactor*DT
        
LPI:    DO I = ML , MU
LPJ:        DO J = NL , NU
    
                h0 = DEP0(I,J)
                K = (h0+TDAMP0)/DZ+OL
                ERO(I,J) = 0.D0 
                tt1 = tt0
                
LPK:            DO
                    taue0 = CRTAU(I,J,K)
                    IF(TDTAU(I,J)>taue0) THEN
                        erosion0 = QE0*(TDTAU(I,J)/taue0-1.D0)*tt1
                        erosion1 = STDEP(K+1)-h0
                        IF(erosion0 <= erosion1) THEN
                            ERO(I,J) = ERO(I,J) + erosion0
                            EXIT LPK
                        ELSE
                            ERO(I,J) = ERO(I,J) + erosion1
                            tt1 = tt1 - erosion1 / QE0 / (TDTAU(I,J)/taue0-1.D0)
                            CRTAU(I,J,K) = CTAUE
                            K = K+1
                            h0 = STDEP(K)
                        ENDIF
                    ELSE
                        EXIT LPK
                    ENDIF

                END DO LPK

            ENDDO LPJ
        ENDDO LPI

#else

        tt0 = MorFactor*DT

LPI:    DO I = ML , MU
LPJ:        DO J = NL , NU
                IF(TDTAU(I,J)>CTAUE) THEN
                    ERO(I,J) = QE0*(TDTAU(I,J)/CTAUE-1.D0)*tt0 !*UPSILON(I,J)
                    TTERO(I,J) = TTERO(I,J) + ERO(I,J)
                ELSE
                    ERO(I,J) = 0.D0
                ENDIF
            ENDDO LPJ
        ENDDO LPI

#endif            

    ! cross-sectional averaged erosion
#ifdef SSC_DIFFUSION_1D
            ERO1D(:) = (2.D0*SUM(ERO,DIM=2) - ERO(:,NL) - ERO(:,NU)) / 2.D0 * DX / BT
#endif

    END SUBROUTINE LocalErosion
 
! - - - - - - - - - - - - - - - - 
!
!       Local Erosion
!
! - - - - - - - - - - - - - - - - 

    SUBROUTINE LocalSSC(Sf,DEP0,SSC1,SSC2,T0,TDEP0)
    
        IMPLICIT NONE
  
        REAL(fp) ,               DIMENSION(:,:) , POINTER , INTENT(INOUT) :: SSC1
        REAL(fp) ,               DIMENSION(:,:) , POINTER , INTENT(INOUT) :: SSC2
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: DEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,           INTENT(IN)    :: Sf
        REAL(fp)                                          , INTENT(IN)    :: T0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(IN)    :: TDEP0
        
        REAL(fp) :: C0 , MSL , DMSL , F1 , F2
        REAL(fp) ,               DIMENSION(:,:) , POINTER                 :: TEMP
        INTEGER :: I , J        

#ifdef SSC_LATERAL_EQ
        DO I = ML , MU
            IF(Sf(I)==0.D0) THEN
                SSC1(I,:) = SSC_CONST
            ELSE
                DO J = NL , NU
                    SSC1(I,J) = SSC_CONST*DEXP(2.D0*WS/SQRT(GRA*Sf(I))/LambdaCZ*(SQRT(1.D0/DEP0(I,NU))-SQRT(1.D0/DEP0(I,J))))
                ENDDO
            ENDIF
        ENDDO

#else
        SSC1(:,:) = SSC_CONST

        !DO I = ML , MU
        !    DO J = NL , NU
        !        IF(TDEP0(I,J)<=0.D0) THEN
        !            SSC2(I,J) = 0.D0
        !        ELSE
        !            C0 = SSC_CONST*DMIN1(2.D0,(TDAMP+DEP0(I,J))) / 2.D0
        !            DMSL = -TDAMP*TDFRQ*DSIN(TDFRQ*T0)
        !            IF(DMSL>0) THEN
        !                F1 = TDEP0(I,J)/DT + DMSL + WS
        !                F2 = TDEP0(I,J)/DT
        !                SSC2(I,J) = (C0*DMSL + F2*SSC1(I,J)) / F1
        !            ELSE
        !                F1 = TDEP0(I,J)/DT + WS
        !                F2 = TDEP0(I,J)/DT
        !                SSC2(I,J) = SSC1(I,J)*F2/F1
        !            ENDIF
        !        ENDIF
        !    END DO
        !END DO
        !TEMP=>SSC2
        !SSC2=>SSC1
        !SSC1=>TEMP       
        
#endif

    END SUBROUTINE LocalSSC
    
! - - - - - - - - - - - -         
!
!   Lateral TRANSPORT
!
! - - - - - - - - - - - -        
    
    SUBROUTINE  LateralST(DEP0,TDTAU,LABL,LAAVA)
    
        IMPLICIT NONE
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: TDTAU
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: DEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: LABL
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: LAAVA
        INTEGER                                                           :: I , J        
        REAL(fp)                                                          :: S0 , S1

        DO I = ML , MU
            DO J = NL , NU           
        
#ifdef BEDLOAD_INC
            IF(TDTAU(I,J)>CTAUE) THEN
                LABL(I,J) = LBFC*(TDTAU(I,J)/CTAUE-1.D0)*DEP0(I,J)**0.3D0 
                !LABL(I,J) = 15.D0*QE0 *(TDTAU(I,J)/CTAUE-1.D0)
            ELSE
                LABL(I,J) = 0.D0
            ENDIF
#endif
            END DO
        END DO
        
        DO I = ML , MU
            DO J = NL , NU-1                
                S0 = (DEP0(I,J+1)-DEP0(I,J))/dx
                S1 = DATAN(S0)

#ifdef AVAlANCHE_INC
                IF(S0>=0.D0) THEN
                    LAAVA(I,J) = DTAN(DMAX1(0.D0,S1-AVARA))*AVAC
                ELSEIF(S0<0.D0) THEN
                    LAAVA(I,J) = DTAN(DMIN1(0.D0,S1+AVARA))*AVAC
                END IF
#endif
            END DO
        END DO
        !CALL Output_Binary_2D('BLB_','.RES',5,LABL)
        !pause
    END SUBROUTINE LateralST

    
    SUBROUTINE AVALANCHE(DEP0)
    
        IMPLICIT NONE
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: DEP0
    
        REAL(fp) , DIMENSION(5000) :: DIFF , AVA
        REAL(fp) :: S0
        INTEGER :: I , J

        DO I = ML , MU
            AVA(:) = 0.D0
            DO J = NL , NU-1
                DIFF(J) = DEP0(I,J+1) - DEP0(I,J)
                S0 = DATAN(DIFF(J)/dx)
                IF(S0>AVARA) THEN
                    AVA(J+1) = AVA(J+1) - 0.4D0*DTAN(AVARA)*dx
                    AVA(J)   = AVA(J)   + 0.4D0*DTAN(AVARA)*dx
                ENDIF
            ENDDO
            DEP0(I,NL:NU) = DEP0(I,NL:NU) + AVA(NL:NU)
        ENDDO
        
    END SUBROUTINE AVALANCHE
! - - - - - - - - - - - -         
!
!   BED EVOLUTION
!
! - - - - - - - - - - - - 
    
! - - - - - 0 order - - - - -- 
    
    SUBROUTINE SSC1D_DIFF_0ORDER(IL,IU,LT,LX,CP1,CP2,GP,DF)
        IMPLICIT NONE
                                                        
        INTEGER  ,                                        INTENT(IN)    :: IL
        INTEGER  ,                                        INTENT(IN)    :: IU
        REAL(fp) ,                                        INTENT(IN)    :: LT
        REAL(fp) ,                                        INTENT(IN)    :: LX
        REAL(fp) ,               DIMENSION(:) , POINTER , INTENT(IN)    :: CP1
        REAL(fp) ,               DIMENSION(:) , POINTER , INTENT(INOUT) :: CP2
        INTEGER  , ALLOCATABLE , DIMENSION(:) ,           INTENT(IN)    :: GP
        REAL(fp) ,                                        INTENT(IN)    :: DF    
        !----------------------------------------------------------------------
        REAL(fp)                                                        :: LXX        
        REAL(fp)                                                        :: DFE , DFW
        INTEGER                                                         :: GDE, GDW
        INTEGER                                                         :: I , IC , K
        
        IC = 0
        LXX = LX**2    
        
        DO I = IL , IU
            
            IC = IC+1
! Local point is dry, set zero
            IF(GP(I)==1) THEN
                amt(1,IC) = 0.D0
                amt(2,IC) = 1.D0
                amt(3,IC) = 0.D0
                amt(4,IC) = 0.D0
                CYCLE  
            ELSEIF(GP(I)==-5) THEN
                amt(1,IC) = 0.D0
                amt(2,IC) = 1.D0
                amt(3,IC) = 0.D0
                amt(4,IC) = CP1(I)
                CYCLE              
! Local point is wet, set diffusion coefficients
            ELSEIF(GP(I)==0) THEN
                !WEST
                IF(I==IL) THEN
                    GDW = 1    ! NEUMANN BOUNDARY
                    !GDW = -2    ! zero Dirichlet
                ELSE
                    GDW = GP(I-1)
                END IF
                IF(GDW==0) THEN
                    DFW = DF * LT / LXX 
                END IF                  
                !EAST                
                IF(I==IU) THEN
                    GDE = 1
                ELSE
                    GDE = GP(I+1)
                END IF                                
                IF(GDE==0) THEN
                    DFE = DF  * LT / LXX
                END IF
!- - - - - - - - - - - - - - - - 
! Left hand side of the equation
!- - - - - - - - - - - - - - - -                 
                IF(GDE==1.AND.GDW<=0) THEN                    
                    amt(1,IC) = - 2.D0*DFW
                    amt(2,IC) = + 2.D0*DFW + 1.D0
                    amt(3,IC) = 0.D0                   
                ELSEIF(GDW==1.AND.GDE<=0) THEN
                    amt(1,IC) =   0.D0
                    amt(2,IC) = + 2.D0*DFE + 1.D0
                    amt(3,IC) = - 2.D0*DFE                    
                ELSEIF(GDW==1.AND.GDE==1) THEN
                    amt(1,IC) =   0.D0
                    amt(2,IC) =   1.D0
                    amt(3,IC) =   0.D0
                ELSEIF(GDW<=0.AND.GDE<=0)  THEN
                    amt(1,IC) = - DFW
                    amt(2,IC) = + DFW + DFE + 1.D0
                    amt(3,IC) = - DFE
                ELSEIF(GDW==-2.AND.GDE<=0)  THEN   ! zero dirichlet boundary
                    amt(1,IC) =   0.D0
                    amt(2,IC) = + DFW + DFE + 1.D0
                    amt(3,IC) = - DFE    
                ELSE
                    WRITE(*,*) I , GP(I-1) , GP(I+1)
                    PAUSE 'tri_matrix_error_1'
                END IF
!- - - - - - - - - - - - - - - - 
! Right hand side of the equation
!- - - - - - - - - - - - - - - - 
                amt(4,IC) = CP1(I)
                
            END IF
            
        END DO

        !if(cd2==15) then
        !OPEN(5,file='amt3.txt')
        !do I = 1, IC
        !    WRITE(5,'(4e25.10,i2)') (amt(K,I),K=1,4),GP(I)
        !end do
        !CLOSE(5)
        !end if
        
        CALL solve3_LU(IC)

        !if(cd2==15) then
        !OPEN(5,file='amt4.txt')
        !do I = 1, IC
        !    WRITE(5,'(10e25.10)') (amt(K,I),K=1,4)
        !end do
        !CLOSE(5)
        !!pause
        !end if
        
        CP2(:) = amt(4,:)

    END SUBROUTINE SSC1D_DIFF_0ORDER  

!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!
!                           UPDATE THE BED ELEVATION
!
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    
    SUBROUTINE BedEvolution_im(TDTAU,CRTAU,STDEP,ERO,SSC1,SSC1D1,LABL,LAAVA,HNDD,EVOSTEP,DEP0)
    
        USE PRISM_RELA    
        IMPLICIT NONE
        
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: TDTAU        
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:,:) ,         INTENT(INOUT) :: CRTAU    
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   ,           INTENT(IN)    :: STDEP  
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: ERO    
        REAL(fp) ,               DIMENSION(:,:) , POINTER , INTENT(INOUT) :: SSC1
        REAL(fp) ,               DIMENSION(:)   , POINTER , INTENT(IN)    :: SSC1D1   
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: LABL           
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: LAAVA
        INTEGER  , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(IN)    :: HNDD
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: EVOSTEP        
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) ,           INTENT(INOUT) :: DEP0 
        !----------------------------------------------------------------------------
        REAL(fp)                                                          :: DEP1 , DEP2
        REAL(fp)                                                          :: CC
        REAL(fp)                                                          :: EROSION0
        REAL(fp)                                                          :: DEPOSITION0
        REAL(fp)                                                          :: BEDLOAD0
        REAL(fp)                                                          :: AVALANCHE0
        INTEGER                                                           :: I , J , JJC , k , K1 , K2
        
! - - - - - - - - - - - - - - - - -
!
!   SUM THE bed evolution
!   Erosion is positive and Deposition is negative
!
! - - - - - - - - - - - - - - - - - 

        DO I = ML , MU
                       
            JJC = 0
            DO J = NL , NU
                JJC = JJC + 1
                
            ! erosion
                    EROSION0 = ERO(I,J)
            ! Deposition
                IF(HNDD(I,J)==0) THEN
#ifdef SSC_DIFFUSION_1D
                        DEPOSITION0 = - SSC1D1(I)*WS*MorFactor*DT
#elif defined(SSC_2D_DIFFUSION)
                        DEPOSITION0 = - SSC1(I)*WS*MorFactor*DT
#elif defined(SSC_CONSTANT)
#ifdef SSC_RELA_PA
                        DEPOSITION0 = - SSC_CONST*WS*PRISMRATIO*MorFactor*DT
#else
                        DEPOSITION0 = - SSC1(I,J)*WS*MorFactor*DT
                        !IF(TDTAU(I,J)>CTAUD) THEN
                        !    DEPOSITION0 = 0 
                        !ELSE
                        !    DEPOSITION0 = -SSC1(I,J)*WS*MorFactor*DT*(-TDTAU(I,J)/CTAUD+1.D0) 
                        !ENDIF
#endif
#endif
                ELSE
                    DEPOSITION0 = 0.D0
                ENDIF

            ! avalanche
#ifdef AVAlANCHE_INC
                IF(J==NL) THEN
                    AVALANCHE0 = LAAVA(I,J)/dx
                ELSEIF(J==NU) THEN
                    AVALANCHE0 = - 2.D0*LAAVA(I,J-1)/dx
                ELSE
                    AVALANCHE0 = (LAAVA(I,J)-LAAVA(I,J-1))/dx
                END IF
#else
                AVALANCHE0 = 0.D0           
#endif
                
                IF(J==NL) THEN
                    AMT(1,JJC) = 0.D0
                    AMT(2,JJC) = 1.D0 + 2.D0*LABL(I,J)*MorFactor*DT/dx**2
                    AMT(3,JJC) =      - 2.D0*LABL(I,J)*MorFactor*DT/dx**2
                    AMT(4,JJC) = DEP0(I,J) + (EROSION0 + DEPOSITION0)              
                ELSEIF(J==NU) THEN
                    AMT(1,JJC) = - 2.D0*LABL(I,J-1)*MorFactor*DT/dx**2
                    AMT(2,JJC) = 1.D0 + 2.D0*LABL(I,J-1)*MorFactor*DT/dx**2
                    AMT(3,JJC) = 0.D0
                    AMT(4,JJC) = DEP0(I,J) + (EROSION0 + DEPOSITION0)
                ELSE
                    AMT(1,JJC) = - LABL(I,J-1)*MorFactor*DT/dx**2
                    AMT(2,JJC) = 1.D0 + (LABL(I,J-1) + LABL(I,J))*MorFactor*DT/dx**2
                    AMT(3,JJC) = - LABL(I,J)*MorFactor*DT/dx**2
                    AMT(4,JJC) = DEP0(I,J) + (EROSION0 + DEPOSITION0)
                ENDIF
           END DO
            
           CALL solve3_LU(JJC)
            
           DO J = 1 , JJC , 1
                DEP1 = DEP0(I,J)
                DEP0(I,J) = AMT(4,J)
                DEP2 = DEP1 + EROSION0
#ifndef TIDE_OSCILLATION
                IF(DEP0(I,J) < 0.005D0) DEP0(I,J)= 0.005D0
#endif                
#ifdef STRATIFICATION
                IF(DEP0(I,J) /= 0.005D0) THEN
                    IF(DEP0(I,J) < DEP1) THEN
                        DO K = OL , OU
                            CRTAU(I,J,K) = DMAX1(CRTAU(I,J,K),CTAUE+(STDEP(K)-DEP0(I,J))*DCRTAU)
                        ENDDO
                    ELSEIF(DEP0(I,J) > DEP1) THEN
                        K1 = INT((DEP2+TDAMP0)/DZ+OL)
                        K2 = INT((DEP0(I,J)+TDAMP0)/DZ+OL)
                        DO K = K1 , K2
                            CRTAU(I,J,K) = DMAX1(CRTAU(I,J,K),CTAUE+(STDEP(K)-DEP0(I,J))*DCRTAU)
                        ENDDO
                    ENDIF
                ENDIF
#endif

           END DO          
                
        END DO    

    END SUBROUTINE BedEvolution_im  
    
    END MODULE Sediment_Transport
    
