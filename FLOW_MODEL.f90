module Flow_Model
    
    use FloatPrecision
    use Parameters
    use Hydro
    use Depth

    implicit none

! Define pointer arrays
    
    REAL(FP) , DIMENSION(:) , POINTER :: QQ1 , QQ2 
    REAL(FP) , DIMENSION(:) , POINTER :: EE1 , EE2
    REAL(FP) , DIMENSION(:) , POINTER :: TEMP    

    contains
    
!- - - - - - - - - - - - - - - - - - -
!
!            Main routine
!
!- - - - - - - - - - - - - - - - - - -
    
    Subroutine Quasi_2D_Flow_Model()
    
end module Flow_Model
    