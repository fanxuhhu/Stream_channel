MODULE LinearSolver
    
    USE FloatPrecision
    
    IMPLICIT NONE
    REAL(fp) , ALLOCATABLE , DIMENSION(:,:) :: AMT
    
    CONTAINS
    
    SUBROUTINE solve5(jc)
        USE FloatPrecision
        IMPLICIT NONE
        INTEGER                      :: jc
        INTEGER                      :: I , J
        REAL(fp)                     :: a

        DO  I=3,jc
            a=amt(3,I-2)
            DO  J=3,6
                amt(J,I-2)=amt(J,I-2)/a
            END DO
            a=-amt(2,I-1)
            amt(3,I-1)=amt(3,I-1)+amt(4,I-2)*a
            amt(4,I-1)=amt(4,I-1)+amt(5,I-2)*a
            amt(6,I-1)=amt(6,I-1)+amt(6,I-2)*a
            a=-amt(1,I)
            amt(2,I)=amt(2,I)+amt(4,I-2)*a
            amt(3,I)=amt(3,I)+amt(5,I-2)*a
            amt(6,I)=amt(6,I)+amt(6,I-2)*a
        END DO
        a=amt(3,jc-1)
        DO  I=3,6
            amt(I,jc-1)=amt(I,jc-1)/a
        END DO
        a=-amt(2,jc)
        amt(3,jc)=amt(3,jc)+a*amt(4,jc-1)
        amt(6,jc)=(amt(6,jc)+a*amt(6,jc-1))/amt(3,jc)
        amt(6,jc-1)=amt(6,jc-1)-amt(4,jc-1)*amt(6,jc)
        DO  I=jc-2,1,-1
            amt(6,I)=amt(6,I)-amt(4,I)*amt(6,I+1)-amt(5,I)*amt(6,I+2)          
        END DO

    END

    SUBROUTINE Solve3_LU(jc)
        
        USE FloatPrecision
        IMPLICIT NONE
        INTEGER                      :: jc
        INTEGER                      :: I , J
        REAL(fp)                     :: temp
        
        DO J = 2 , jc
            temp = amt(1,J) / amt(2,J-1)
            amt(1,J) = temp
            amt(2,J) = amt(2,J) - temp * amt(3,J-1)
        END DO
        DO J = 2 , jc
            amt(4,J) = amt(4,J) - amt(4,J-1) * amt(1,J)
        END DO
        amt(4,jc) = amt(4,jc) / amt(2,jc)
        DO J = jc-1 , 1 , -1
            amt(4,J) = (amt(4,J) - amt(3,J)*amt(4,J+1)) / amt(2,J)
        END DO

    END SUBROUTINE solve3_LU

    SUBROUTINE SOLVE3(JC)
        USE FloatPrecision
        IMPLICIT NONE
        INTEGER                      :: JC
        INTEGER                      :: I , J
        REAL(fp)                     :: A

        DO  I=2,JC
            A=AMT(2,I-1)
            DO  J=2,4
                AMT(J,I-1)=AMT(J,I-1)/A
            END DO
            A=-AMT(1,I)
            AMT(2,I)=AMT(2,I)+AMT(3,I-1)*A
            AMT(4,I)=AMT(4,I)+AMT(4,I-1)*A
        END DO
        A=AMT(2,JC)
        DO  J=2,4
           AMT(J,JC)=AMT(J,JC)/A
        END DO

        A=AMT(3,JC-1)
        DO  I=JC-1,1,-1
            AMT(4,I)=AMT(4,I)-AMT(3,I)*AMT(4,I+1)
        END DO
        
        RETURN
    END SUBROUTINE SOLVE3
    
END MODULE LinearSolver