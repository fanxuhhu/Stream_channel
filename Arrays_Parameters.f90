#include "AMYCAI.h"

!= = = = = = = = = = = = = = = = = = = = = = = = = = 
!
! modules in this file defines all the global Parameters
! AND arrays in the model...
!
!= = = = = = = = = = = = = = = = = = = = = = = = = = 

!- - - - - - - - - - - - 
!    double precision
!- - - - - - - - - - - -
    
    MODULE FloatPrecision

        IMPLICIT NONE
    
        INTEGER, parameter :: fp = selected_real_kind(8)
    
    END MODULE FloatPrecision

!- - - - - - - - - - - - 
!    double precision
!- - - - - - - - - - - -    

    MODULE TimeCatcher

        USE FloatPrecision
    
        IMPLICIT NONE
        
        REAL(fp):: T0
        INTEGER :: CD , CD2 , CD3
        INTEGER :: IDT , IT0 , IHOUR , IDAY , IYEAR
    
    END MODULE TimeCatcher

!- - - - - - - - - - - - 
!    up AND low bounds
!    of arrays
!- - - - - - - - - - - -

    MODULE ArraySize

        IMPLICIT NONE
    
        INTEGER , TARGET :: ML , MU , NL , NU , OL , OU
    
    END MODULE ArraySize

    MODULE HYDRO
        USE FloatPrecision
        IMPLICIT NONE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET :: Q1   , Q2
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET :: U1   , U2
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET :: ETA1 , ETA2 , DETA1

    END MODULE HYDRO
    
    MODULE Topography
        USE FloatPrecision
        IMPLICIT NONE
        
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: MAREA
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)   :: TDUD   , TDTAU , RHOGSD
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET :: FC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: Sf
        REAL(fp) , ALLOCATABLE , DIMENSION(:)            :: EtaU        
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: XDistance
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: Dep0_GOHST
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)  :: TDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: DEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:,:) :: CRTAU
        REAL(fp) , ALLOCATABLE , DIMENSION(:) :: STDEP
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: CRTAUL
        INTEGER  , ALLOCATABLE , DIMENSION(:)   :: piecewise_right
        INTEGER  , ALLOCATABLE , DIMENSION(:)   :: piecewise_right_Q
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: Dep1
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: DepQ
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: DDEP
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: DDepQ
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) , TARGET :: Upsilon
        
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: BED_CURVATURE
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: SUBMERGE
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: BED_CURVATURE_Q
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: EVO
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: AC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: BC
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: HR
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: BCEVO
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: MDEPEVO
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: MDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: MDEPU
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: MDEPE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: MTDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: MTDEPU
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: MTDEPE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   :: HMTDEP0
        INTEGER  , ALLOCATABLE , DIMENSION(:)   :: ND
        INTEGER  , ALLOCATABLE , DIMENSION(:)   :: NDQ
        INTEGER  , ALLOCATABLE , DIMENSION(:)   :: WetSegments
        REAL(fp)                                :: DCRI

        REAL(fp) , ALLOCATABLE , DIMENSION(:)    :: Q_RESIDUAL1       
        REAL(fp) , ALLOCATABLE , DIMENSION(:)    :: Q_RESIDUAL2          
        REAL(fp) , ALLOCATABLE , DIMENSION(:)    :: HB        
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)  :: HTDEP0
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)  :: SSC_DIFFU_HWL_X
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)  :: SSC_DIFFU_HWL_Y
        REAL(fp) , ALLOCATABLE , DIMENSION(:)    :: HDIFF
        INTEGER  , ALLOCATABLE , DIMENSION(:)    :: HND
        INTEGER  , ALLOCATABLE , DIMENSION(:,:)  :: HNDD
        
        REAL(fp) , ALLOCATABLE , DIMENSION(:)    :: HETA
                
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           :: MaxTau , MeanTau
        
    END MODULE Topography

    MODULE Suspended_Load
        USE FloatPrecision
        IMPLICIT NONE
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET :: SC1D1, SC1D2
        REAL(fp) , ALLOCATABLE , DIMENSION(:)            :: ERO1D
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) , TARGET :: SC2D1, SC2D2
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:) , TARGET :: EROQ, ERO
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)          :: TTERO
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)          :: EVOSTEP
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)          :: Evolution_TDCYC      
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)          :: LABL
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)          :: LAAVA
        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)          :: Lateral_Transport_Q
    End MODULE Suspended_Load
    
    MODULE ReGrid
        USE FloatPrecision
        IMPLICIT NONE
        REAL(fp)                                         :: Re_dx
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET :: Re_TTDep
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET :: Re_Upsilon0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET :: Re_U0
        REAL(fp) , ALLOCATABLE , DIMENSION(:)   , TARGET :: Re_TAU0
        INTEGER  , ALLOCATABLE , DIMENSION(:)            :: JSTART
        INTEGER  , ALLOCATABLE , DIMENSION(:)            :: JEND
    END MODULE ReGrid
    
!- - - - - - - - - - - - - - - 
!   spatial-temporal step size 
!- - - - - - - - - - - - - - - 

    MODULE GRIDS

        USE FloatPrecision
    
        IMPLICIT NONE
    
        REAL(fp) , TARGET :: DX , DY , DZ , DT
        LOGICAL           :: BND_SPEC
    
    END MODULE GRIDS

!- - - - - - - - - - - - - - - 
!    Global parameterws
!- - - - - - - - - - - - - - - 
    
    MODULE Parameters

        USE FloatPrecision
    
        IMPLICIT NONE

        REAL(fp)             :: BT , LU , HV
        REAL(fp)             :: Pcs
        REAL(fp)             :: Pbl
        REAL(fp)             :: Ini_Longi_slope
        REAL(fp)             :: Ini_land_Dep
        REAL(fp)             :: Pertub
        REAL(fp)             :: AVARA
        REAL(fp)             :: AVAC
        REAL(fp)             :: MorFactor
        REAL(fp)             :: TD_Excursion_X
        REAL(fp)             :: TD_Excursion_Y
        REAL(fp)             :: Dasterisk
        REAL(fp)             :: DCRTAU
        REAL(fp)             :: CTAUE
        REAL(fp)             :: CTAUD
        REAL(fp)             :: CTAUB
        REAL(fp)             :: QE0
        REAL(fp)             :: SSC_CONST
        REAL(fp)             :: WS
        REAL(fp)             :: DDRY
        REAL(fp)             :: DWET
        REAL(fp)             :: Lambda
        REAL(fp)             :: LambdaCZ
        REAL(fp)             :: TDAMP0
        REAL(fp)             :: TDAMP
        REAL(fp)             :: MANNING
        REAL(fp)             :: EDDYVITY
        REAL(fp)             :: D50
        REAL(fp)             :: Ki_Visc
        REAL(fp)             :: LBFC
        REAL(fp)             :: rhoRsgd
        REAL(fp)             :: Q_CONST
        INTEGER              :: YEARS
        INTEGER              :: HistoryID
        INTEGER              :: Start_File_ID
        INTEGER              :: End_File_ID
        INTEGER              :: File_ID_FRQ
        INTEGER              :: File_Number
        REAL(fp) , PARAMETER :: Dsand   = 62.D-6
        REAL(fp) , PARAMETER :: GRA     = 9.8D0
        REAL(fp) , PARAMETER :: PI      = 4.D0*DATAN(1.D0)
        REAL(fp) , PARAMETER :: RHOW    = 1000.D0
        REAL(fp) , PARAMETER :: RHOS    = 2600.D0
        REAL(fp) , PARAMETER :: REDY    = (RHOS-RHOW)/RHOW
        REAL(fp) , PARAMETER :: TDCYC   = 43200.D0
        REAL(fp) , PARAMETER :: TDFRQ   = 2.D0*PI/TDCYC
        INTEGER  , PARAMETER :: DAYS    = 360 
        INTEGER  , PARAMETER :: HOURS   = 24
        INTEGER  , PARAMETER :: SECONDS = 3600
    
    END MODULE Parameters
    
    
    MODULE PRISM_RELA
        USE FloatPrecision
        IMPLICIT NONE    
        REAL(fp) :: TDPRISM , TDPRISM0 , TDPRISM1, PRISMRATIO, TDAREA
    END MODULE