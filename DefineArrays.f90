#include "AMYCAI.h"

SUBROUTINE Define_Arrays
USE ArraySize
USE Topography
USE HYDRO
USE Suspended_Load
USE ReGrid
USE LinearSolver

IMPLICIT NONE
INTEGER :: status
CHARACTER(len=80) :: err_msg

! 1D HYDRO ARRAIES READY
ALLOCATE(Q1(ML:MU), stat=status,errmsg=err_msg)
ALLOCATE(Q2(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(U1(ML:MU), stat=status,errmsg=err_msg)
ALLOCATE(U2(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(EtaU(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(ETA1(ML:MU+1),stat=status,errmsg=err_msg)
ALLOCATE(DETA1(ML:MU+1),stat=status,errmsg=err_msg)
ALLOCATE(ETA2(ML:MU+1),stat=status,errmsg=err_msg)
ALLOCATE(BC(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(AC(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(HR(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(BCEVO(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(MDEPEVO(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(HB(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(Sf(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(Q_RESIDUAL1(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(Q_RESIDUAL2(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(FC(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(MDEP0(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(MTDEP0(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(HMTDEP0(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(MDEPU(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(MTDEPU(ML:MU),  stat=status,errmsg=err_msg)
ALLOCATE(MDEPE(ML:MU+1),  stat=status,errmsg=err_msg)
ALLOCATE(MTDEPE(ML:MU+1),  stat=status,errmsg=err_msg)
ALLOCATE(MAREA(ML:MU),  stat=status,errmsg=err_msg)


ALLOCATE(Re_TTDep(2000),  stat=status,errmsg=err_msg)
ALLOCATE(Re_Upsilon0(2000),  stat=status,errmsg=err_msg)
ALLOCATE(Re_U0(2000),  stat=status,errmsg=err_msg)
ALLOCATE(Re_TAU0(2000),  stat=status,errmsg=err_msg)
ALLOCATE(JSTART(2000),  stat=status,errmsg=err_msg)
ALLOCATE(JEND(2000),  stat=status,errmsg=err_msg)


! 2D HYDRO ARRAIES READY
ALLOCATE(TDUD(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(TDTAU(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(RHOGSD(ML:MU,NL:NU),stat=status,errmsg=err_msg)
! DepTH ARRAIES READ

ALLOCATE(WetSegments(ML:MU),stat=status,errmsg=err_msg)
ALLOCATE(MaxTau(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(MeanTau(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(XDistance(NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(DEP0(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(CRTAU(ML:MU,NL:NU,OL:OU),stat=status,errmsg=err_msg)
ALLOCATE(CRTAUL(NL:NU,OL:OU),stat=status,errmsg=err_msg)
ALLOCATE(STDEP(OL:OU),stat=status,errmsg=err_msg)
ALLOCATE(TDEP0(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(DEP0_GOHST(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(BED_CURVATURE(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(SUBMERGE(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(BED_CURVATURE_Q(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(piecewise_right(ML:MU),stat=status,errmsg=err_msg)
ALLOCATE(piecewise_right_Q(ML:MU),stat=status,errmsg=err_msg)
ALLOCATE(Dep1(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(HTDEP0(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(SSC_DIFFU_HWL_X(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(HDIFF(ML:MU),stat=status,errmsg=err_msg)
ALLOCATE(SSC_DIFFU_HWL_Y(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(DDEP(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(LABL(ML:MU,NL:NU),stat=status,errmsg=err_msg)
LABL(:,:) = 0.D0
ALLOCATE(LAAVA(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(Lateral_Transport_Q(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(DepQ(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(DDepQ(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(Upsilon(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(EVO(ML:MU,NL:NU) ,stat=status,errmsg=err_msg)
ALLOCATE(HND(ML:MU),stat=status,errmsg=err_msg)
ALLOCATE(HNDD(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(ND(ML:MU+1),stat=status,errmsg=err_msg)
ALLOCATE(NDQ(ML:MU),stat=status,errmsg=err_msg)
ALLOCATE(HETA(ML:MU),stat=status,errmsg=err_msg)


! Suspended Load Arrays Ready
ALLOCATE(ERO1D(ML:MU),stat=status,errmsg=err_msg)
ALLOCATE(ERO(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(TTERO(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(EROQ(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(SC2D1(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(SC2D2(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(SC1D1(ML:MU),stat=status,errmsg=err_msg)
ALLOCATE(SC1D2(ML:MU),stat=status,errmsg=err_msg)
ALLOCATE(EVOSTEP(ML:MU,NL:NU),stat=status,errmsg=err_msg)
ALLOCATE(Evolution_TDCYC(ML:MU,NL:NU),stat=status,errmsg=err_msg)

! LinearSolver
ALLOCATE(AMT(6,2*(MAX0(MU,NU)+100)),stat=status,errmsg=err_msg)

END SUBROUTINE Define_Arrays