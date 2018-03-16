!Morphodynamic routine or calculate the PA relationship
#define MORPHODYNAMIC_ROUTINE
!#define PA_ROUTINE

!Calculate long-term evolution or short-term flow, SSC field 
#define LONGTERM
#ifndef LONGTERM
#define SHORTTERM
#endif 

!READ history Dep file
!#define HISTORYREAD

!Flow Model
!#define SAINT_VENENT
#ifndef SAINT_VENENT
#define CONTINUITY
#endif

!tidal amplitude
#define CONST_AMP
!#define INCREASE_AMP
!#define DYNAMIC_AMP

!If tide oscillation include (only for continuity flow model)
#ifdef CONTINUITY
!#define TIDE_OSCILLATION
!#define TIDE_MAX_DISCHARGE
#endif

!firiction formula
#define Shiono_Knight_1988

!Wave
!#define WAVEINC

!choose sediment transport model

! artificial evolution
!#define ARTIFICIAL_EVO

!1D
!#define SSC_DIFFUSION_1D
!#define RESIDUALINC

! or 2D
!#define SSC_DIFFUSION_2D

! or constant
#define SSC_CONSTANT
!#define SSC_LATERAL_EQ
!#define SSC_RELA_PA

! cut off SSC
!#define CUT_OFF_SSC 2900

! Lateral transport
!#define AVALANCHE2
!#define BEDLOAD_INC
!#define ARTIFICIAL_LBFC

!
!#define MAXTAUEVO

! transport modes
#define DALPAOS

! sediment stratification
!#define STRATIFICATION

! wave
!#define WAVE_INC
