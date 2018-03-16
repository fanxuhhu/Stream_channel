#include "AMYCAI.h"
    
SUBROUTINE BATHYMETRY
    
    USE FloatPrecision
    USE Parameters
    USE ArraySize
    USE Topography
    USE GRIDS
    IMPLICIT NONE
    INTEGER :: I , J , K
    INTEGER :: status , ISTAT
    CHARACTER(len=80) :: err_msg
    CHARACTER(len=14) :: HistoryFile
    REAL(fp) :: y, y1, y2, z1, z2, slope1, slope2, amp0 , PTBR
    INTEGER , DIMENSION(1) :: SEED = (/3/)
    REAL(fp) :: RANDOMDep
    INTEGER :: JJ

! - - - - - - - - - - - - - - - - 
!
! Calculate the width and length of the 
! computational basin. We calculate half
! of the channel, because the channel is
! symmetry! 
!
! - - - - - - - - - - - - - - - - 
    
    BT = DBLE(NU-NL)*dx
    LU = DBLE(MU-ML)*dy
    
! - - - - - - - - - - - - - - - - 
!
! Define the initial shape of the
! computational basin
!
! - - - - - - - - - - - - - - - -
    
! = = = = = = = = = = = = = = 
!  Linear slope case
! = = = = = = = = = = = = = =

    DO I = MU , ML , -1
        y = DBLE(MU-I)*dy
        DEP0(I,:) = Ini_land_Dep + y*Ini_Longi_slope
    END DO
    
    amp0 = Pertub
    slope2 = amp0/BT ! Half of the basin
    
    DO J = NL , NU , 1
        !IF(MOD(J,2)==0) THEN
        !    PTBR = 0.1
        !ELSE
        !    PTBR = -0.1
        !ENDIF    
        DEP0(:,J) = DEP0(:,J) + slope2*DBLE(J-NL)*dx !+ PTBR
    END DO  
    
! - - - - - - - - - - -    
!   READ FROM HISTORY
! - - - - - - - - - - -

#ifdef HISTORYREAD

    WRITE(HistoryFile,'(A4,I6.6,A4)') 'Dep_',HistoryID,'.RES'
    OPEN(unit=10,file=HistoryFile,form='unformatted',status='old',iostat=istat,action='READ')
    READ(10) DEP0
    CLOSE(10)
    
#endif


#ifdef STRATIFICATION

    DO K = OL , OU
        STDEP(K) = - TDAMP0 + DZ*DBLE(K-OL)
    ENDDO
    
    DO I = ML , MU
        DO J = NL , NU
            DO K = OL , OU
                CRTAU(I,J,K) = DMAX1(0.4D0,0.4D0+(STDEP(K)-DEP0(I,J))*DCRTAU)
            ENDDO
        ENDDO
    ENDDO
    
#endif

! - - - - - - - - - - -    
!   Output initial topography
! - - - - - - - - - - -

    OPEN(UNIT=10,FILE='BATHYMETRY0.TXT' ,form='unformatted',status='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
    WRITE(10) DEP0
    CLOSE(10)

! - - - - - - - - - - -    
! Calculate Lateral distance
! - - - - - - - - - - -   
        
        DO J = NL , NU
            XDistance(J) = dble(J-1)*dx
        END DO

    END SUBROUTINE BATHYMETRY
    
    
! = = = = = = = = = = = = = = 
!  discontinuous Linear slope case
! = = = = = = = = = = = = = =
    
    ! Define the longitudinal shape
    !
    !                       Land region z2 = 0.5m
    !                      ------------
    !                     /
    !    z1=50m Dep      /  Intertidal region 
    !                   /   linear slope 0.1%
    !       Outer Sea  /
    !     -------------   
    !
    !    0-----------y1---y2--------->LU
    !
    != = = = = = = = = = = = = = = = 
    !
    !    y1 = 20000.D0
    !    y2 = 25000.D0    
    !
    !DO I = ML , MU
    !    y = DBLE(I-1)*dy
    !    IF(y <= y1) THEN
    !        DEP0(I,:) = z1
    !    ELSEIF(y>y1.AND.y<y2) THEN
    !        DEP0(I,:) = z1 - (y-y1)*slope1
    !    ELSEIF(y>=y2) THEN
    !        DEP0(I,:) = z2 !- (y-y2)*slope2
    !    ENDIF
    !END DO

!  cosine incision
   
    !JJ = 1
    !DO J  = NL , NU , 1
    !    IF(J>=JJ.AND.J<=NU+1-JJ) THEN
    !        DEP0(:,J) = DEP0(:,J) - amp0*COS(2.D0*PI/DBLE(NU+1-2*JJ)*DBLE(J-JJ)) + amp0
    !    END IF
    !END DO

! only erosion case

    !DEP0(:,:) = -1.5D0*TDAMP0
    !z1 = 1.5D0*TDAMP0
    !z2 = -TDAMP0
    !slope1 = (z1-z2)/DBLE(MU-ML)/dy
    !JJ = 3
    !DO I = ML , MU
    !    y = DBLE(I-ML)*dy
    !    DEP0(I,(NL+NU)/2-JJ:(NL+NU)/2+JJ) = z1 - slope1*y
    !END DO
    
! = = = = = = = = = = = = = = 
!  van de wegen slope case
! = = = = = = = = = = = = = =

    !DO I = MU , ML , -1
    !    y = DBLE(MU-I)*dy
    !    IF(y<=5000.D0) THEN
    !        DEP0(I,:) = Ini_land_Dep + y*Ini_Longi_slope
    !    ELSE
    !        DEP0(I,:) = 1000.D0
    !    END IF
    !END DO
    !
    !amp0 = 0.01D0
    !
    !
    !slope2 = 2.D0*amp0/DBLE(NU-NL)*dx
    !
    !DO J = NL , (NL+NU)/2 , 1
    !    DEP0(:,J) = DEP0(:,J) + slope2*DBLE(J-NL)*dx
    !    DEP0(:,NL+NU-J) = DEP0(:,J)
    !END DO    
