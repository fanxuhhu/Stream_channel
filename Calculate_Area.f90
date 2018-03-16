!Subroutine Calculate_Area
!    use FloatPrecision
!    use ArraySize
!    use grids
!    use depth
!    use hydro
!    implicit none
!
!    integer :: i , j
!    real(fp), dimension(NL:NU) :: darea
!
!    do i = ML , MU
!        do j = NL , NU
!            darea(j) = dmax1(0.d0,ttdep(i,j))
!        end do
!        A1(i) = (2.d0*sum(darea(:))-darea(NL)-darea(NU))*dx/2.d0
!    end do
!
!End Subroutine Calculate_Area