!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
!!
!! This program is desined to simulate the morphodynamics of a OPEN estuary or tidal channel.
!! It is not a triditional 2D model that solve fullformed shallow water equations OR  
!! sediment transport equation. The Flow model is composed by two mutually coupled 1D
!! models: a Longitudinal 1D model which solves Saint Venent euations AND a transversal 
!! 1D model which solves the transversal distribution of shear stress....
!!
!! This program is writen by Xu Fan in 2016, when he was studying in the University of
!! Auckland, NZ, with the guidance of Giovanni Coco.
!! 
!! This program is writen for Cai Duo ....
!!    
!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
    
#include "AMYCAI.h"

program AMYCAI
    
    IMPLICIT NONE
    
    CALL ReadData    
    CALL Define_Arrays   

! define the sediment transport mode
#ifdef DALPAOS
    CALL DALPAOS_SETTINGS
#elif defined(DAVIES)
    CALL DAVIES_SETTINGS
#elif defined(IKEDA)
    CALL IKEDA_SETTINGS
#endif
    
#ifdef MORPHODYNAMIC_ROUTINE

    CALL Morphodynamic
    
#elif defined(PA_ROUTINE)

    CALL PA_RELATION

#endif

END program AMYCAI

