! = = = = = = = = = = = = = = = = 
!
!  Diffusion proccess for SSC
!
! = = = = = = = = = = = = = = = =     
   
!    SUBROUTINE SSC_Diffusion_2D_Model(SSC1,SSC2)
!
!        IMPLICIT NONE
!
!        REAL(fp) , DIMENSION(:,:) , POINTER :: SSC1 , SSC2 , TEMP
!        INTEGER :: I , J , i1 , i2 , j1 , j2
!
!
!! Calculate local erosion
!
!        CALL LocalErosion
!
!!- - - - - - - - - - - - - - 
!!    ADI method is used
!!- - - - - - - - - - - - - - 
!    
!! Y-Direcsion
! 
!        LX   => DY
!        LY   => DX
!        LT   => DT
!        LXX = LX**2
!        LYY = LY**2
!        
!        ILOW => ML
!        IUP  => MU
!        JLOW => NL
!        JUP  => NU
!        
!        DO J = JLOW , JUP
!            
!            j1 = J+1
!            j2 = J-1
!            IF(J==JLOW) j2=J
!            IF(J==JUP) j1=J
!            
!! results of last step
!            
!            VARS  => SSC1(:,j2)
!            VARP  => SSC1(:,J)
!            VARN  => SSC1(:,j1)
!            
!! wet AND dry information
!            
!            GIDS  => HNDD(:,j2)
!            GIDP  => HNDD(:,J)
!            GIDN  => HNDD(:,j1)
!            
!! diffusion coefficient -Y
!            
!            VSCP  => SSC_DIFFU_HWL_Y(:,J)
!            
!! diffusion coefficient -X
!            
!            VSCS  => SSC_DIFFU_HWL_X(:,j2)
!            VSCN  => SSC_DIFFU_HWL_X(:,J)
!            
!! erosion term: local source term
!            
!            CEP   => ERO(:,J)
!            
!! local water Depth
!            
!            HP    => HTDEP0(:,J)
!
!! result
!            
!            FRUIT => SSC2(:,J)     
!
!! ADI solver
!            
!            CALL Tridiagonal_ADI('Y',J)
!            
!        END DO
!
!! Exchange pointers    
!    
!        TEMP=>SSC2
!        SSC2=>SSC1
!        SSC1=>TEMP
!        NULLIFY(VARS,VARP,VARN,FRUIT,GIDS,GIDP,GIDN,VSCP,VSCS,VSCN,TEMP,CEP) 
!        NULLIFY(LX,LY,LT,ILOW,IUP,JLOW,JUP)
! 
! ! - - - - - - - - - - - - - - - - -
!
! ! X-Direcsion
!
!        LX   => DX
!        LY   => DY
!        LT   => DT
!        LXX = LX**2
!        LYY = LY**2
!        
!        ILOW => NL
!        IUP  => NU
!        JLOW => ML
!        JUP  => MU
!        
!        DO J = JLOW , JUP
!            
!            j1 = J+1
!            j2 = J-1
!            IF(J==JLOW) j2=J
!            IF(J==JUP ) j1=J
!            
!! results of last step
!            
!            VARS  => SSC1(j2,:)
!            VARP  => SSC1(J,:)
!            VARN  => SSC1(j1,:)
!            
!! wet AND dry information
!            
!            GIDS  => HNDD(j2,:)
!            GIDP  => HNDD(J,:)
!            GIDN  => HNDD(j1,:)
!            
!! diffusion coefficient -Y
!            
!            VSCP  => SSC_DIFFU_HWL_X(J,:)
!            
!! diffusion coefficient -X
!            
!            VSCS  => SSC_DIFFU_HWL_Y(j2,:)
!            VSCN  => SSC_DIFFU_HWL_Y(J,:)
!            
!! erosion term: local source term
!            
!            CEP   => ERO(J,:)
!            
!! local water Depth
!            
!            HP    => HTDEP0(J,:)
!
!! result
!            
!            FRUIT => SSC2(J,:)     
!
!! ADI solver
!            
!            CALL Tridiagonal_ADI('X',J)
!            
!        END DO
!
!! Exchange pointers    
!    
!        TEMP=>SSC2
!        SSC2=>SSC1
!        SSC1=>TEMP
!        NULLIFY(VARS,VARP,VARN,FRUIT,GIDS,GIDP,GIDN,VSCP,VSCS,VSCN,TEMP,CEP) 
!        NULLIFY(LX,LY,LT,ILOW,IUP,JLOW,JUP)
!
!    END SUBROUTINE SSC_Diffusion_2D_Model
    
! = = = = = = = = = = = = = = = =    
!
!     SSC constant model
!
! = = = = = = = = = = = = = = = =  
    
    !SUBROUTINE SSC_Constant_Model(SSC1)
    !    IMPLICIT NONE
    !
    !    REAL(fp) , DIMENSION(:,:) , POINTER :: SSC1
    !    
    !    CALL LocalErosion
    !
    !END SUBROUTINE SSC_Constant_Model
    
    
    ! 2D routine
    
!    SUBROUTINE Bed_Evolution(SSC1)
!    
!        IMPLICIT NONE
!        
!        REAL(fp) , DIMENSION(:,:) , POINTER :: SSC1
!        INTEGER :: I , J
!        
!! Erosion is positive and Deposition is negative
!        
!        DO I = ML , MU
!            DO J = NL , NU
!                EVOSTEP(I,J) = ERO(I,J) - SSC1(I,J)*WS
!                Evolution_TDCYC(I,J) = Evolution_TDCYC(I,J) + EVOSTEP(I,J)
!            END DO
!        END DO
!        
!    END SUBROUTINE Bed_Evolution

    
    
! = = = = = = = = = = = = = = = = = = = = 
!
!   ADI ROUTINE: solve tridiagonal matrix
!
! = = = = = = = = = = = = = = = = = = = =     

!    SUBROUTINE Tridiagonal_ADI(DIR,ROW)    
!    
!        IMPLICIT NONE
!        
!        CHARACTER(LEN=1) :: DIR
!        INTEGER          :: ROW
!
!        REAL(fp) :: AMT(6,3001)
!        REAL(fp) :: DCE , DCW , DCN , DCS
!        INTEGER  :: I , IC , K
!        INTEGER  :: IDE, IDW, IDS, IDN
!        IC = 0
!        
!        DO I = ILOW , IUP
!            
!            IC = IC+1
!
!! Local point is dry
!
!            IF(GIDP(I)==1) THEN
!
!! set zero
!                
!                amt(1,IC) = 0.D0
!                amt(2,IC) = 1.D0
!                amt(3,IC) = 0.D0
!                amt(4,IC) = 0.D0
!                CYCLE
!                
!! Local point is wet
!                
!            ELSEIF(GIDP(I)==0) THEN
!
!!- - - - - - - - - - - - - -
!! Set diffusion coefficients             
!!- - - - - - - - - - - - - -
!                
!!SOUTH
!                
!                IF(ROW==JLOW) THEN
!                    IDS = 1
!                ELSE
!                    IDS = GIDS(I)
!                END IF
!
!                IF(IDS==0) THEN
!                    DCS = VSCS(I)  * LT / LYY * 0.5D0 / HP(I)
!                END IF  
!
!!NORTH
!                
!                IF(ROW==JUP) THEN
!                    IDN = 1
!                ELSE
!                    IDN = GIDN(I)
!                END IF
!
!                IF(IDN==0) THEN
!                    DCN = VSCN(I)  * LT / LYY * 0.5D0 / HP(I)
!                END IF  
!
!!WEST
!
!                IF(I==ILOW) THEN
!                    IDW = 1
!                ELSE
!                    IDW = GIDP(I-1)
!                END IF
!
!                IF(IDW==0) THEN
!                    DCW = VSCP(I-1) * LT / LXX * 0.5D0 / HP(I)  
!                END IF  
!                
!!EAST                
!
!                IF(I==IUP) THEN
!                    IDE = 1
!                ELSE
!                    IDE = GIDP(I+1)
!                END IF                
!                
!                IF(IDE==0) THEN
!                    DCE = VSCP(I)  * LT / LXX * 0.5D0 / HP(I)
!                END IF
!                
!!- - - - - - - - - - - - - - - - 
!! Left hand side of the equation
!!- - - - - - - - - - - - - - - -                 
!
!                IF(IDE==1.AND.IDW==0) THEN
!                    
!                    amt(1,IC) = - 2.D0*DCW
!                    amt(2,IC) = + 2.D0*DCW + 1.D0
!                    amt(3,IC) = 0.D0
!                    
!                ELSEIF(IDW==1.AND.IDE==0) THEN
!                    
!                    amt(1,IC) =   0.D0
!                    amt(2,IC) = + 2.D0*DCE + 1.D0
!                    amt(3,IC) = - 2.D0*DCE
!                    
!                ELSEIF(IDW==1.AND.IDE==1) THEN
!                    
!                    amt(1,IC) =   0.D0
!                    amt(2,IC) =   1.D0
!                    amt(3,IC) =   0.D0
!                    
!                ELSEIF(IDW==0.AND.IDE==0)  THEN
!                    
!                    amt(1,IC) = - DCW
!                    amt(2,IC) = + DCW + DCE + 1.D0
!                    amt(3,IC) = - DCE
!                    
!                ELSE
!                    
!                    WRITE(*,*) ROW, GIDP(I-1) , GIDP(I+1)
!                    PAUSE 'tri_matrix_error_1'
!                    
!                END IF
!                
!!- - - - - - - - - - - - - - - - 
!! Left hand side of the equation
!!- - - - - - - - - - - - - - - - 
!                
!                IF(    IDN==1.AND.IDS==0) THEN
!                    
!                    amt(4,IC) = VARP(I) + 2.D0*DCS*(VARS(I)-VARP(I))
!                    
!                ELSEIF(IDS==1.AND.IDN==0) THEN
!                    
!                    amt(4,IC) = VARP(I) + 2.D0*DCN*(VARN(I)-VARP(I))
!                    
!                ELSEIF(IDS==1.AND.IDN==1) THEN
!                    
!                    amt(4,IC) = VARP(I)
!                    
!                ELSEIF(IDS==0.AND.IDN==0)  THEN
!                    
!                    amt(4,IC) = VARP(I) + DCN*(VARN(I)-VARP(I)) + DCS*(VARS(I)-VARP(I))
!
!                ELSE
!                    
!                    WRITE(*,*) ROW, IDS , IDN
!                    PAUSE 'tri_matrix_error_2'
!                    
!                END IF
!                
!        ! sourse term
!                
!                amt(4,IC) = amt(4,IC) + CEP(I)/HP(I)*0.5D0*LT
!                amt(2,IC) = amt(2,IC) + WS/HP(I)*0.5D0*LT
!                
!            END IF
!            
!        END DO
!
!    !    OPEN(5,file='amt.txt')
!    !    do I = 1, IC
!    !        WRITE(5,'(10e25.10)') (amt(K,I),K=1,4)
!    !    end do
!    !    CLOSE(5)
!
!        CALL solve3_LU(IC)
!
!        fruit(:) = amt(4,:)
!        
!        !IF(pros=='sedi') THEN
!        !IF(dir=='y') THEN
!        !OPEN(5,file='amt2.txt')
!        !do I = 1, IC
!        !    WRITE(5,'(10e25.10)') (amt(K,I),K=1,4)
!        !end do
!        !CLOSE(5)
!        !pause 'xg'
!        !end IF
!
!    END SUBROUTINE Tridiagonal_ADI     

    
!- - - - - - - - - - - - - -
!    
!   Bed evolution per step
!
!- - - - - - - - - - - - - -

! 1D routine
