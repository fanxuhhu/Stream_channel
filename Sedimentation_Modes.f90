#include "AMYCAI.h"
    
SUBROUTINE DALPAOS_SETTINGS
    USE FloatPrecision
    USE Grids
    USE Topography
    USE Parameters
    USE ArraySize
    IMPLICIT NONE
    
! volumetric
    SSC_CONST = SSC_CONST / RHOS
    QE0 = QE0 / RHOS
    
! LBFC  Lateral Bedload flux coefficient
#ifndef ARTIFICIAL_LBFC

    IF(D50<Dsand) THEN
        LBFC = 1.D0/18.D0*Dsand*D50**(-0.1D0)*(REDY*GRA)**(0.4D0)*Ki_Visc**(0.2D0) * Pbl
    ELSE
        LBFC = 1.D0/18.D0*D50  *D50**(-0.1D0)*(REDY*GRA)**(0.4D0)*Ki_Visc**(0.2D0) * Pbl
    ENDIF

#endif
! CTAUB
    Dasterisk = (REDY*GRA/Ki_Visc**2)**(1.D0/3.D0)*D50
    IF(D50 < Dsand) THEN
        ! van Rijn 2007
        IF(Dasterisk < 4.D0) THEN
            CTAUB = 0.115D0*Dasterisk**(-0.5D0)*(RHOS-RHOW)*GRA*D50*(Dsand/D50)
        ELSE
            CTAUB = 0.14D0*Dasterisk**(-0.64D0)*(RHOS-RHOW)*GRA*D50*(Dsand/D50)
        ENDIF
    ELSE
        IF(Dasterisk < 4.D0) THEN
            CTAUB = 0.115D0*Dasterisk**(-0.5D0)*(RHOS-RHOW)*GRA*D50*(1.D0+Pcs)**3
        ELSE
            CTAUB = 0.14D0*Dasterisk**(-0.64D0)*(RHOS-RHOW)*GRA*D50*(1.D0+Pcs)**3
        ENDIF
    ENDIF
    ! soulsby 1997
    !CTAUB = GRA*(RHOS-RHOW)*D50*(0.3D0/(1.D0+1.2D0*Dasterisk)+0.055D0*(1.D0-dexp(-0.02D0*Dasterisk)))

    CTAUB = CTAUE

! ctaub
    
    
! radian
    AVARA = AVARA * PI / 180.D0    
    
END SUBROUTINE DALPAOS_SETTINGS
    
SUBROUTINE IKEDA_SETTINGS
    USE FloatPrecision
    USE Grids
    USE Topography
    USE Parameters
    USE ArraySize
    IMPLICIT NONE
    
! volumetric
    SSC_CONST = SSC_CONST / RHOS

! LBFC  Lateral Bedload flux coefficient
    LBFC = 1.D0/18.D0*Dsand*D50**(-0.1D0)*(REDY*GRA)**(0.4D0)*Ki_Visc**(0.2D0)
    
! radian
    AVARA = AVARA * PI / 180.D0  
    
! rhoRsgd
    rhoRsgd = RHOW*REDY*GRA*D50

! QE0
    QE0 = 0.001D0 / D50 / GRA / WS / RHOW**2 / REDY*0.001D0

! LBFC  Lateral Bedload flux coefficient
    LBFC = 1.85D0 / D50**(1.5D0) / GRA**(1.5D0) / RHOW**2 / REDY**(1.5D0)*0.001D0

END SUBROUTINE IKEDA_SETTINGS    
    