#include "AMYCAI.h"
    
SUBROUTINE ReadData
    USE FloatPrecision
    USE Grids
    USE Topography
    USE Parameters
    USE ArraySize
    IMPLICIT NONE
    INTEGER  :: ISTAT , I

    OPEN(UNIT=10,FILE='input_amycai.dat',status='OLD',IOSTAT=ISTAT,ACTION='READ')

!CHECK OPEN ERROR

    IN_OK: IF(ISTAT/=0) THEN
        WRITE(*,'(A41)')    '|---------------------------------------'
        WRITE(*,'(A35,I3)') '|  Input file OPEN error, IOSTAT =' , ISTAT
        WRITE(*,'(A41)')    '|  Please check the <<input_amycai.dat>>'
        WRITE(*,'(A39)')    '|  Report from SUBROUTINE <<ReadData>>'
        WRITE(*,'(A41)')    '|---------------------------------------'
        PAUSE
    END IF IN_OK

! READ DATAS

    READ(10,*) ML
    READ(10,*) MU
    READ(10,*) NL
    READ(10,*) NU
    READ(10,*) OL
    READ(10,*) OU
    READ(10,*) DX
    READ(10,*) DY
    READ(10,*) DZ
    READ(10,*) DT
    READ(10,*) TDAMP0
    READ(10,*) Q_CONST
    READ(10,*) YEARS
    READ(10,*) MANNING
    READ(10,*) lambda
    READ(10,*) LambdaCZ
    READ(10,*) EDDYVITY
    READ(10,*) DDRY
    READ(10,*) DWET
    READ(10,*) TD_Excursion_X
    READ(10,*) TD_Excursion_Y
    READ(10,*) CTAUE
    READ(10,*) CTAUD
    READ(10,*) DCRTAU
    READ(10,*) Pcs
    READ(10,*) Pbl
    READ(10,*) LBFC
    READ(10,*) QE0
    READ(10,*) SSC_CONST
    READ(10,*) WS
    READ(10,*) DCRI
    READ(10,*) Ki_Visc
    READ(10,*) D50
    READ(10,*) Ini_Longi_slope
    READ(10,*) Ini_land_Dep
    READ(10,*) Pertub
    READ(10,*) MorFactor
    READ(10,*) AVARA
    READ(10,*) AVAC
    READ(10,*) HistoryID
    READ(10,*) Start_File_ID
    READ(10,*) End_File_ID
    READ(10,*) File_ID_FRQ
    CLOSE(10)

!! LBFC  Lateral Bedload flux coefficient
!    LBFC = 1.D0/18.D0*Dsand*D50**(-0.1D0)*(REDY*GRA)**(0.4D0)*Ki_Visc**(0.2D0)
!! radian
!    AVARA = AVARA * PI / 180.D0    
!    
!! CTAUB
!    Dasterisk = (REDY*GRA/Ki_Visc**2)**(1.D0/3.D0)*D50
!    !CTAUE = GRA*(RHOS-RHOW)*D50*(0.3D0/(1.D0+1.2D0*Dasterisk)+0.055D0*(1.D0-dexp(-0.02D0*Dasterisk)))
!
!! WS
!    !WS = Ki_Visc/D50*(SQRT(10.36D0**2+1.049D0*Dasterisk**3)-10.36D0)
!    !write(*,*) WS
!    !pause
!
!! volumetric SSC_CONST
!    SSC_CONST = SSC_CONST / RHOS
    
    
! FILE_NUMBER : PA module, continuity flow model
    File_Number = (End_File_ID-Start_File_ID)/File_ID_FRQ+1
    
END SUBROUTINE ReadData