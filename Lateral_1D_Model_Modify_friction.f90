#include "AMYCAI.h"

MODULE Lateral_1D_Model

    USE TimeCatcher
    USE FloatPrecision
    USE ArraySize
    USE GRIDS
    USE hydro
    USE Topography
    USE Parameters

    IMPLICIT NONE

    contains

SUBROUTINE Cross_Sectional_Tau(Surface,Discharge)

    IMPLICIT NONE
    
    REAL(fp) , ALLOCATABLE , POINTER , DIMENSION(:)   :: Surface
    REAL(fp) , ALLOCATABLE , POINTER , DIMENSION(:)   :: Discharge
    
    REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           :: TTDepthQ
    REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           :: DWFriction 
    REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           :: LateralDiff
    REAL(fp) , ALLOCATABLE , DIMENSION(:)             :: Sf
    REAL(fp) , DIMENSION(6,3001)                      :: amt
    
    INTEGER                                           :: Regrid_Time , NG
    REAL(fp) , parameter                              :: tol = 1.d-6     
    REAL(fp)                                          :: CC1, CC2, fff, dc 
    INTEGER                                           :: i, j, j1, j2, k, jc 
    INTEGER                                           :: status
    CHARACTER(len=80)                                 :: err_msg

!   regrid--
    
    REAL(fp) , ALLOCATABLE , DIMENSION(:)           :: Re_Ud        
    REAL(fp) , ALLOCATABLE , DIMENSION(:)           :: Re_TAU_TD        
    REAL(fp) , ALLOCATABLE , DIMENSION(:)           :: Re_Upsilon    
    REAL(fp) , ALLOCATABLE , DIMENSION(:)           :: Re_TTDepthQ
    REAL(fp) , ALLOCATABLE , DIMENSION(:)           :: Re_DWFriction 
    REAL(fp) , ALLOCATABLE , DIMENSION(:)           :: Re_LateralDiff    
    REAL(fp)                                        :: Re_dx , Re_x , xl , xr
    INTEGER                                         :: jj
    INTEGER  , ALLOCATABLE , DIMENSION(:)           :: JSTART , JEND
 
!   effective width
    INTEGER  , ALLOCATABLE , DIMENSION(:)           :: GAMA_JSTART , GAMA_JEND
    REAL(fp) , ALLOCATABLE , DIMENSION(:)           :: ttdep
    REAL(fp)                                        :: cc
    REAL(fp)                                        :: aa1 , aa2
    REAL(fp)                                        :: bb1 , bb2 , bbb
    INTEGER                                         :: piecewise0 , piecewise1
    
! ALLOCATE arrays

    ALLOCATE(ttdep(NL:NU),stat=status,errmsg=err_msg)
    ALLOCATE(TTDepthQ(ML:MU,NL:NU),stat=status,errmsg=err_msg)
    ALLOCATE(DWFriction(ML:MU,NL:NU),stat=status,errmsg=err_msg)
    ALLOCATE(LateralDiff(ML:MU,NL:NU),stat=status,errmsg=err_msg)
    ALLOCATE(Sf(ML:MU),stat=status,errmsg=err_msg)
    ALLOCATE(JSTART(ML:MU),stat=status,errmsg=err_msg)
    ALLOCATE(JEND(ML:MU),stat=status,errmsg=err_msg)
    ALLOCATE(GAMA_JSTART(ML:MU),stat=status,errmsg=err_msg)
    ALLOCATE(GAMA_JEND(ML:MU),stat=status,errmsg=err_msg)    

! - - - - - - - - - - - - - - - - 
! get the water depth of the cross-section
! get the lateral slope cross-section
! - - - - - - - - - - - - - - - - 
    
    DO i = ML , MU-1

! upwind scheme
        
        IF(Discharge(i)>=0.D0) THEN
            TTDepthQ(i,:) = DEP0(i,:) + Surface(i)            
        ELSE IF(Discharge(i)<0.D0) THEN
            TTDepthQ(i,:) = DEP0(i+1,:) + Surface(i+1)
        END IF
        
! slope and upsilon
        
        DO j = NL , NU
            IF(j==NL) THEN
                DDEPQ(i,j) = (TTDepthQ(i,j+1)-TTDepthQ(i,j)) / dx
            ELSEIF(j==NU) THEN
                DDEPQ(i,j) = (TTDepthQ(i,j)-TTDepthQ(i,j-1)) / dx
            ELSE
                DDEPQ(i,j) = (TTDepthQ(i,j+1)-TTDepthQ(i,j-1)) / 2.D0 / dx
            END IF
            UPSILON(i,j) = sqrt(1+DDEPQ(i,j)**2)
        END DO        
        
    END DO
    
! - - - - - - - - - - - - - - - -     
! get the wet points of each cross-section
! - - - - - - - - - - - - - - - - 
    
    DO i = ML , MU-1
        JSTART(i) = NL
        JEND(i) = NU
        WetPoints(i) = 0
        LPJ9: DO j = NL , NU
            IF(TTDepthQ(i,j)>0.D0) THEN
                JSTART(i) = j-1
                JEND(i)  = NL + NU - JSTART(i)
                WetPoints(i) = -JSTART(i) + JEND(i) - 1
                EXIT LPJ9
            END IF          
        END DO LPJ9
        IF(JSTART(i)>(NU+NL)/2) THEN
            PAUSE 'error regrid'
        END IF
    END DO

!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
!
!   LPI: Iterate the cross-sectional shear stress...
!
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    
    LPI:DO i = ML , MU-1


        
        TAU_TD(i,:) = 0.D0
        Ud(i,:)  = 0.D0
        gama(i) = 0.D0       
        
! When discharge is zero, shear stress AND velocity are both zero. 
        
        IF(Discharge(i)>0.D0)  THEN
            piecewise0 = NL+NU-Piecewise_Right(i)
            piecewise1 = Piecewise_Right(i)
        ELSE IF(Discharge(i)<0.D0) THEN
            piecewise0 = NL+NU-Piecewise_Right(i+1)
            piecewise1 = Piecewise_Right(i+1)
        ELSE IF(Discharge(i)==0.D0)  THEN
            cycle LPI            
        END IF        
        
! = = = = = = = = = = 
!
! To determine if re-griding is needed       
!
! = = = = = = = = = = 
        
        IF (WetPoints(i)>=19) THEN
        
        
! Calculate the friction coefficients
        
            DO j = NL , NU
                IF(TTDepthQ(i,j)>0.D0) THEN
                    DWFriction(i,j) = manning**2*8.D0*gra*TTDepthQ(i,j)**(-1.D0/3.D0)
                ELSE
                    DWFriction(i,j) = 0.D0
                END IF
            END DO  
        
! Calculate the diffusion term
        
            DO j = NL , NU-1
                IF(TTDepthQ(i,j)>0.D0.AND.TTDepthQ(i,j+1)>0.D0) THEN
                    fff = manning**2*8.D0*gra*(0.5d0*(TTDepthQ(i,j)+TTDepthQ(i,j+1)))**(-1.D0/3.D0)
                    LateralDiff(i,j) = 0.125D0*Lambda*(TTDepthQ(i,j)+TTDepthQ(i,j+1))**2*sqrt(fff/8.D0)
                ELSE
                    LateralDiff(i,j) = 0.D0
                END IF
            END DO
        
! initial guess of Sf : Surface slope
        
            Sf(i) = 0.5d0
        
!  - - - - - - - - - - - - - - - - - - - - - -
! LPK: the iteration loop
! The upper limit of iteration times is 1000
!- - - - - - - - - - - - - - - - - - - - - - -
        
            LPK: DO k = 1 , 1000
 
                jc = 0

! - - - - - - - - - - - - - - - - - - - - - - - - 
! LPJ: the loop calculate the transversal distribution 
!      of the shear stress AND the velocity
! - - - - - - - - - - - - - - - - - - - - - - - -        
    
                LPJ:DO j = NL , NU
                
                    j1 = j+1
                    j2 = j-1
                    IF(j==NL) j2 = j
                    IF(j==NU) j1 = j
                    jc = jc + 1

! Dry point
                
                    IF(TTDepthQ(i,j)<=0) THEN
                        amt(1,jc) = 0.D0
                        amt(2,jc) = 1.D0
                        amt(3,jc) = 0.D0
                        amt(4,jc) = 0.D0
                        cycle LPJ
                    END IF

! Edge points

                    IF(j==NL.OR.j==NU) THEN
                        amt(1,jc) = 0.D0       
                        amt(2,jc) = DWFriction(i,j)/8.D0*UPSILON(i,j) 
                        amt(3,jc) = 0.D0 
                        amt(4,jc) = gra*TTDepthQ(i,j)*dabs(Sf(i))
                        cycle LPJ
                    END IF

! General case, first calculate the diffusion terms
                
                    CC1 = 1.D0/dx**2*LateralDiff(i,j-1)     
                    CC2 = 1.D0/dx**2*LateralDiff(i,j)

! Fill the band matrix

                    amt(1,jc) = - CC2        
                    amt(2,jc) = DWFriction(i,j)/8.D0*UPSILON(i,j) + CC1 + CC2 
                    amt(3,jc) = - CC1 
                    amt(4,jc) = gra*TTDepthQ(i,j)*dabs(Sf(i))

                END DO LPJ

!            IF(iiday==3099) THEN
                !open(111,file='amt1.txt')
                !DO j = 1 , jc
                !    write(111,'(50e25.15)') (amt(jj,j),jj=1,4),DWFriction(i,j)/8.D0,LateralDiff(I,J),UPSILON(i,j)
                !END DO
                !close(111)

!            END IF
                CALL solve3_LU(jc,amt)
            
    !        IF(i==57.AND.iicdxf==1099) THEN
                !open(111,file='amt2.txt')
                !DO j = 1 , jc
                !    write(111,'(50e25.15)') (amt(jj,j),jj=1,4)
                !END DO
                !close(111)
                !pause 'amt6'
    !        END IF
            
! get the velocity

                DO j = NL , NU
                    Ud(i,j) = sqrt(amt(4,j))
                    IF(Discharge(i)<0.D0) Ud(i,j) = - Ud(i,j)
                END DO

! estimate the discharge

                dc = 0.D0
                DO j = NL , NU-1
                !DO j = piecewise0 , piecewise1-1
                    dc = dc + 0.25D0*(Ud(i,j)+Ud(i,j+1))*(dmax1(0.D0,TTDepthQ(i,j+1))+dmax1(0.D0,TTDepthQ(i,j)))*dx
                END DO

! compare with the actual discharge

                IF(dabs(dc-Discharge(i)) <= tol) THEN
                    DO j = 1 , jc
                        TAU_TD(i,j) = DWFriction(i,j)/8.D0*Ud(i,j)**2*rhow  !*UPSILON(i,j)
                    END DO
                    
! calculate gama
! Coupling the lateral 1d model with the longitudinal model.
! i.e., modify the friction coefficient of the longitudinal
! 1d model to match the intergral bed resistance calculated
! by 1d lateral model.
                    
                    gama(i) = 0.D0

                    DO j = NL , NU-1 , 1
                        IF(DWFriction(i,j)>0.D0.OR.DWFriction(i,j+1)>0.D0) THEN
                            gama(i) = gama(i) &
                                    & + 0.5D0 * (&
                                    &   DWFriction(i,j)/8.D0*Ud(i,j)**2*UPSILON(i,j) &
                                    & + DWFriction(i,j+1)/8.D0*Ud(i,j+1)**2*UPSILON(i,j+1) &
                                    & ) * dx
                        END IF
                    END DO
                    gama(i) = gama(i) / dabs(Discharge(i))
    
                    cycle LPI
                END IF
                                
                Sf(i) = Sf(i)*(Discharge(i)/dc)**2

            END DO LPK

! error message

            IF(k==1001) THEN
                write(*,*) i , CD
                pause 'error1'
            END IF
            
! = = = = = = = = 
! = = = = = = = = 
!            
! IF Wet Points is too few ----- Regriding
!            
! = = = = = = = =
! = = = = = = = = 
            
        ELSEIF(WetPoints(i)<19) THEN

! calculate multiple times
            
            Regrid_Time = 21 / (WetPoints(i)+2) + 1

! number of new points

            NG = Regrid_Time*WetPoints(i) + Regrid_Time + 1

! New Grid size
            
            Re_dx = dx / DBLE(Regrid_Time) 
            
! Allocate arrays
            
            ALLOCATE(Re_TTDepthQ(1:NG),stat=status,errmsg=err_msg)
            ALLOCATE(Re_DWFriction(1:NG),stat=status,errmsg=err_msg)
            ALLOCATE(Re_LateralDiff(1:NG),stat=status,errmsg=err_msg)
            ALLOCATE(Re_Upsilon(1:NG),stat=status,errmsg=err_msg)
            ALLOCATE(Re_Ud(1:NG),stat=status,errmsg=err_msg)
            ALLOCATE(Re_TAU_TD(1:NG),stat=status,errmsg=err_msg)
            
! Calculate TTDepthQ after regrid
    
            DO j = 1 , NG
                Re_x = xx(JSTART(i)) + dble(j-1)*Re_dx      
                LPR1:DO jj = JSTART(i) , JEND(i)-1
                    xl = xx(jj)
                    xr = xx(jj+1)
                    IF(Re_x>=xl.and.Re_x<=xr) THEN
                        !IF(jj==JSTART(i)) THEN ! extrapolate
                            xl = xx(jj+1)
                            xr = xx(jj+2)                            
                            Re_TTDepthQ(j) = ((Re_x-xl)*TTDepthQ(i,jj+2)+(xr-Re_x)*TTDepthQ(i,jj+1))/dx
                            Re_Upsilon(j) = ((Re_x-xl)*Upsilon(i,jj+2)+(xr-Re_x)*Upsilon(i,jj+1))/dx
                        !ELSEIF(jj==JEND(i)-1) THEN ! extrapolate
                        !    xl = xx(jj-1)
                        !    xr = xx(jj)
                        !    Re_TTDepthQ(j) = ((Re_x-xl)*TTDepthQ(i,jj)+(xr-Re_x)*TTDepthQ(i,jj-1))/dx
                        !    Re_Upsilon(j) = ((Re_x-xl)*Upsilon(i,jj)+(xr-Re_x)*Upsilon(i,jj-1))/dx
                        !ELSE ! inpolate
                        !    Re_TTDepthQ(j) = ((Re_x-xl)*TTDepthQ(i,jj+1)+(xr-Re_x)*TTDepthQ(i,jj))/dx
                        !    Re_Upsilon(j) = ((Re_x-xl)*Upsilon(i,jj+1)+(xr-Re_x)*Upsilon(i,jj))/dx
                        !ENDIF
                        
                        EXIT LPR1
                    END IF
                END DO LPR1
            END DO
            !Re_TTDepthQ(1) = TTDepthQ(i,JSTART(i)) 
            !Re_TTDepthQ(NG) =TTDepthQ(i,JEND(i)) 
            
! Calculate the friction coefficients after regrid
        
            DO j = 1 , NG
                IF(Re_TTDepthQ(j)>0.D0) THEN
                    Re_DWFriction(j) = manning**2*8.D0*gra*Re_TTDepthQ(j)**(-1.D0/3.D0)
                ELSE
                    Re_DWFriction(j) = 0.D0
                END IF
            END DO  
        
! Calculate the diffusion term after regrid
        
            DO j = 1 , NG-1
                IF(Re_TTDepthQ(j)>0.D0.AND.Re_TTDepthQ(j+1)>0.D0) THEN
                    fff = manning**2*8.D0*gra*(0.5d0*(Re_TTDepthQ(j)+Re_TTDepthQ(j+1)))**(-1.D0/3.D0)
                    Re_LateralDiff(j) = 0.125D0*Lambda*(Re_TTDepthQ(j)+Re_TTDepthQ(j+1))**2*sqrt(fff/8.D0)
                ELSE
                    Re_LateralDiff(j) = 0.D0
                END IF
            END DO

! initial guess of Sf : Surface slope
        
            Sf(i) = 0.5d0
        
!  - - - - - - - - - - - - - - - - - - - - - -
! LPK: the iteration loop
! The upper limit of iteration times is 1000
!- - - - - - - - - - - - - - - - - - - - - - -
        
            LPK2: DO k = 1 , 1000
 
                jc = 0

! - - - - - - - - - - - - - - - - - - - - - - - - 
! LPJ: the loop calculate the transversal distribution 
!      of the shear stress AND the velocity
! - - - - - - - - - - - - - - - - - - - - - - - -        
    
                LPJ2:DO j = 1 , NG
                
                    j1 = j+1
                    j2 = j-1
                    IF(j==1) j2 = j
                    IF(j==NG) j1 = j
                    jc = jc + 1

! Dry point
                
                    IF(Re_TTDepthQ(j)<=0) THEN
                        amt(1,jc) = 0.D0
                        amt(2,jc) = 1.D0
                        amt(3,jc) = 0.D0
                        amt(4,jc) = 0.D0
                        cycle LPJ2
                    END IF

! Edge points

                    !IF(j==1.OR.j==NG) THEN
                    !    write(*,*) cd,i
                    !    write(*,*) wetpoints(i),NG
                    !    write(*,*) jstart(i),jend(i)
                    !    write(*,*) Re_TTDepthQ(1),Re_TTDepthQ(NG)
                    !    write(*,*) TTDepthQ(i,jstart(i)),TTDepthQ(i,jend(i))                        
                    !    pause 'regrid error2'
                    !END IF

! General case, first calculate the diffusion terms
                
                    CC1 = 1.D0/Re_dx**2*Re_LateralDiff(j-1)     
                    CC2 = 1.D0/Re_dx**2*Re_LateralDiff(j)

! Fill the band matrix

                    amt(1,jc) = - CC2        
                    amt(2,jc) = Re_DWFriction(j)/8.D0*Re_UPSILON(j) + CC1 + CC2 
                    amt(3,jc) = - CC1 
                    amt(4,jc) = gra*Re_TTDepthQ(j)*dabs(Sf(i))

                END DO LPJ2

!            IF(iiday==3099) THEN
                !open(111,file='amt1.txt')
                !DO j = 1 , jc
                !    write(111,'(50e25.15)') (amt(jj,j),jj=1,4),DWFriction(i,j)/8.D0,LateralDiff(I,J),UPSILON(i,j)
                !END DO
                !close(111)

!            END IF
                CALL solve3_LU(jc,amt)
            
    !        IF(i==57.AND.iicdxf==1099) THEN
                !open(111,file='amt2.txt')
                !DO j = 1 , jc
                !    write(111,'(50e25.15)') (amt(jj,j),jj=1,4)
                !END DO
                !close(111)
                !pause 'amt6'
    !        END IF
            
! get the velocity

                DO j = 1 , NG
                    Re_Ud(j) = sqrt(amt(4,j))
                    IF(Discharge(i)<0.D0) Re_Ud(j) = - Re_Ud(j)
                END DO

! estimate the discharge

                dc = 0.D0
                DO j = 1 , NG-1
                    dc = dc + 0.25D0*(Re_Ud(j)+Re_Ud(j+1))*(dmax1(0.D0,Re_TTDepthQ(j+1))+dmax1(0.D0,Re_TTDepthQ(j)))*Re_dx
                END DO

! compare with the actual discharge

                IF(dabs(dc-Discharge(i)) <= tol) THEN
                    
                    DO j = 1 , jc
                        Re_TAU_TD(j) = Re_DWFriction(j)/8.D0*Re_Ud(j)**2*rhow
                    END DO
                    
! fit back to coarse grid
                    
                    jj = JSTART(i) - 1
                    DO j = 1 , NG , Regrid_Time
                        jj = jj + 1
                        Ud(i,jj) = Re_Ud(j)
                        TAU_TD(i,jj) = Re_TAU_TD(j)
                    END DO
 
! calculate gama
                    
                    gama(i) = 0.D0
                    DO j = 1 , NG-1
                        IF(Re_DWFriction(j)>0.D0.OR.Re_DWFriction(j+1)>0.D0) THEN
                            gama(i) = gama(i) &
                                    & + 0.5D0 * (&
                                    &   Re_DWFriction(j)/8.D0*Re_Ud(j)**2*Re_UPSILON(j) &
                                    & + Re_DWFriction(j+1)/8.D0*Re_Ud(j+1)**2*Re_UPSILON(j+1) &
                                    & ) * Re_dx
                        END IF
                    END DO
                    gama(i) = gama(i) / dabs(Discharge(i))
                    
!DEALLOCATE
                    
                    DEALLOCATE(Re_TTDepthQ,stat=status,errmsg=err_msg)
                    DEALLOCATE(Re_DWFriction,stat=status,errmsg=err_msg)
                    DEALLOCATE(Re_LateralDiff,stat=status,errmsg=err_msg)
                    DEALLOCATE(Re_Upsilon,stat=status,errmsg=err_msg)
                    DEALLOCATE(Re_Ud,stat=status,errmsg=err_msg)
                    DEALLOCATE(Re_TAU_TD,stat=status,errmsg=err_msg)   
                    
                    cycle LPI
                    
                END IF
                
                Sf(i) = Sf(i)*(Discharge(i)/dc)**2

            END DO LPK2

! error message

            IF(k==1001) THEN
                write(*,*) i , CD , wetpoints(i)
                pause 'LPK error1'
            END IF
            
        END IF
            
    END DO LPI


! DEALLOCATE arrays

    DEALLOCATE(ttdep,stat=status,errmsg=err_msg)
    DEALLOCATE(TTDepthQ,stat=status,errmsg=err_msg)
    DEALLOCATE(DWFriction,stat=status,errmsg=err_msg)
    DEALLOCATE(LateralDiff,stat=status,errmsg=err_msg)
    DEALLOCATE(Sf,stat=status,errmsg=err_msg)
    DEALLOCATE(JSTART,stat=status,errmsg=err_msg)
    DEALLOCATE(JEND,stat=status,errmsg=err_msg)   
    DEALLOCATE(GAMA_JSTART,stat=status,errmsg=err_msg)
    DEALLOCATE(GAMA_JEND,stat=status,errmsg=err_msg)   
    
END SUBROUTINE Cross_Sectional_Tau

END MODULE Lateral_1D_Model