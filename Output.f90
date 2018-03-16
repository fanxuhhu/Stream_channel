#include "AMYCAI.h"
    
    MODULE OUTPUT
    
        USE ArraySize
        USE FloatPrecision
        IMPLICIT NONE
                
    CONTAINS
    
!- - - - - - 
!
! Binary_2D
!
!- - - - - - 

    SUBROUTINE Output_Binary_3D(Head,Tail,ID,RES)
        
        IMPLICIT NONE

        REAL(fp) , ALLOCATABLE , DIMENSION(:,:,:)           , INTENT(IN)   :: RES
        CHARACTER(LEN=4)                  :: HEAD , TAIL
        INTEGER                           :: ID


        CHARACTER(LEN=14) :: FILENAME
        INTEGER           :: ISTAT
    
        WRITE(FILENAME,'(A4,I6.6,A4)') HEAD,ID,TAIL
        OPEN(UNIT=10,FILE=FILENAME,FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
        WRITE(10) RES
        CLOSE(10)

    END SUBROUTINE OUTPUT_BINARY_3D
    
    SUBROUTINE Output_Binary_2D(Head,Tail,ID,RES)
        
        IMPLICIT NONE

        REAL(fp) , ALLOCATABLE , DIMENSION(:,:)           , INTENT(IN)   :: RES
        CHARACTER(LEN=4)                  :: HEAD , TAIL
        INTEGER                           :: ID


        CHARACTER(LEN=14) :: FILENAME
        INTEGER           :: ISTAT
    
        WRITE(FILENAME,'(A4,I6.6,A4)') HEAD,ID,TAIL
        OPEN(UNIT=10,FILE=FILENAME,FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
        WRITE(10) RES
        CLOSE(10)

    END SUBROUTINE OUTPUT_BINARY_2D

    SUBROUTINE Output_Binary_2D_P(Head,Tail,ID,RES)
        
        IMPLICIT NONE

        REAL(fp) ,               DIMENSION(:,:)   , POINTER , INTENT(INOUT) :: RES
        CHARACTER(LEN=4)                  :: HEAD , TAIL
        INTEGER                           :: ID


        CHARACTER(LEN=14) :: FILENAME
        INTEGER           :: ISTAT
    
        WRITE(FILENAME,'(A4,I6.6,A4)') HEAD,ID,TAIL
        OPEN(UNIT=10,FILE=FILENAME,FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
        WRITE(10) RES
        CLOSE(10)

    END SUBROUTINE OUTPUT_BINARY_2D_P
    
    
 !- - - - - - 
!
! Binary_1D
!
!- - - - - - 
 
    SUBROUTINE Output_Binary_1D(Head,Tail,ID,RES) 
        
        IMPLICIT NONE

        REAL(fp) , ALLOCATABLE , DIMENSION(:)           , INTENT(IN)   :: RES
        CHARACTER(LEN=4)                  :: HEAD , TAIL
        INTEGER                           :: ID


        CHARACTER(LEN=14) :: FILENAME
        INTEGER           :: ISTAT
    
        WRITE(FILENAME,'(A4,I6.6,A4)') HEAD,ID,TAIL
        OPEN(UNIT=10,FILE=FILENAME,FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
        WRITE(10) RES
        CLOSE(10)
    
    END SUBROUTINE OUTPUT_BINARY_1D

!- - - - - - 
!
! Binary_1D  INTEGER
!
!- - - - - - 

    SUBROUTINE Output_Binary_1D_INTEGER(Head,Tail,ID,RES)

        IMPLICIT NONE

        INTEGER , ALLOCATABLE , DIMENSION(:)           , INTENT(IN)   :: RES
        CHARACTER(LEN=4)                  :: HEAD , TAIL
        INTEGER                           :: ID


        CHARACTER(LEN=14) :: FILENAME
        INTEGER           :: ISTAT
    
        WRITE(FILENAME,'(A4,I6.6,A4)') HEAD,ID,TAIL
        OPEN(UNIT=10,FILE=FILENAME,FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
        WRITE(10) RES
        CLOSE(10)
    
    END SUBROUTINE OUTPUT_BINARY_1D_INTEGER

    
END MODULE OUTPUT
!
!SUBROUTINE OUTPUT(HEAD,TAIL,ID,RES,I1,I2,IP,J1,J2,JP)
!    USE ArraySize
!    USE FloatPrecision
!    CHARACTER*4 :: HEAD , TAIL
!    INTEGER     :: ID
!    REAL(fp) , DIMENSION(ML:MU,NL:NU) :: RES
!    INTEGER :: I1 , I2 , IP , J1 , J2 , JP
!    INTEGER :: I , J
!    CHARACTER*14 :: FILENAME
!    INTEGER :: status , ISTAT
!    CHARACTER(len=80) :: err_msg
!    
!    WRITE(FILENAME,'(A4,I6.6,A4)') HEAD,ID,TAIL
!    OPEN(UNIT=10,FILE=FILENAME,status='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
!    DO I = I1 , I2 , IP
!        WRITE(10,'(2000E25.10)') (RES(I,J),J=J1,J2,JP)
!    END DO
!    CLOSE(10)
!    
!END SUBROUTINE OUTPUT
!
!SUBROUTINE OUTPUTINT(HEAD,TAIL,ID,RES,I1,I2,IP,J1,J2,JP)
!    USE ArraySize
!    USE FloatPrecision
!    CHARACTER*4 :: HEAD , TAIL
!    INTEGER     :: ID
!    INTEGER , DIMENSION(ML:MU,NL:NU) :: RES
!    INTEGER :: I1 , I2 , IP , J1 , J2 , JP
!    INTEGER :: I , J
!    CHARACTER*14 :: FILENAME
!    INTEGER :: status , ISTAT
!    CHARACTER(len=80) :: err_msg
!    
!    WRITE(FILENAME,'(A4,I6.6,A4)') HEAD,ID,TAIL
!    OPEN(UNIT=10,FILE=FILENAME,status='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
!    DO I = I1 , I2 , IP
!        WRITE(10,'(2000I3)') (RES(I,J),J=J1,J2,JP)
!    END DO
!    CLOSE(10)
!    
!    END SUBROUTINE OUTPUTINT
!    
!    
!SUBROUTINE OUTPUT1D(HEAD,TAIL,ID,RES,I1,I2,IP)
!    USE ArraySize
!    USE FloatPrecision
!    CHARACTER*4 :: HEAD , TAIL
!    INTEGER     :: ID
!    REAL(fp) , DIMENSION(ML:MU) :: RES
!    INTEGER :: I1 , I2 , IP
!    INTEGER :: I , J
!    CHARACTER*14 :: FILENAME
!    INTEGER :: status , ISTAT
!    CHARACTER(len=80) :: err_msg
!    
!    WRITE(FILENAME,'(A4,I6.6,A4)') HEAD,ID,TAIL
!    OPEN(UNIT=10,FILE=FILENAME,status='REPLACE',IOSTAT=ISTAT,ACTION='WRITE')
!    DO I = I1 , I2 , IP
!        WRITE(10,'(200E25.10)') RES(I)
!    END DO
!    CLOSE(10)
!    
!    END SUBROUTINE OUTPUT1d
!    
    
