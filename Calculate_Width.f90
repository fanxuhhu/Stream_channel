SUBROUTINE Calculate_Width
    USE FloatPrecision
    USE ArraySize
    USE Parameters
    USE grids
    USE hydro
    USE Topography
    IMPLICIT NONE

    INTEGER :: i , j
    REAL(fp), DIMENSION(NL:NU) :: darea
    REAL(fp) :: c1 , c2 , maxdd

    c1 = 1.3D0
    c2 = 0.9D0
    DO i = ML , MU
        maxdd = c1*TDAMP + dmin1(maxval(DEP0(i,:)),c2*TDAMP)
        IF(maxdd > 0.D0) THEN
            DO j = NL , NU
                darea(j) = dmax1(0.D0,(c1*TDAMP + dmin1(DEP0(i,j),c2*TDAMP)))
            END DO
            B1(i) = (2.D0*sum(darea(:))-darea(NL)-darea(NU))*dx/2.D0/maxdd
        ELSE
            B1(i) = 0.D0 !
        END IF
    END DO

END SUBROUTINE Calculate_Width