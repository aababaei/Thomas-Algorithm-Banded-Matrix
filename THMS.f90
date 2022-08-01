      PROGRAM EXAMPLE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(7,8)

!     Example system:
!     |-6  1  1  1  0  0  0 |  3  |
!     | 1 -6  1  1  1  0  0 |  1  |
!     | 1  1 -6  1  1  1  0 |  0  |
!     | 1  1  1 -6  1  1  1 |  0  |
!     | 0  1  1  1 -6  1  1 | -8  |
!     | 0  0  1  1  1 -6  1 | -17 |
!     | 0  0  0  1  1  1 -6 | -27 |

!     Band rearrangement:
      T(1,:) = (/ 0, 0, 0, -6, 1, 1, 1, 3 /)
      T(2,:) = (/ 0, 0, 1, -6, 1, 1, 1, 1 /)
      T(3,:) = (/ 0, 1, 1, -6, 1, 1, 1, 0 /)
      T(4,:) = (/ 1, 1, 1, -6, 1, 1, 1, 0 /)
      T(5,:) = (/ 1, 1, 1, -6, 1, 1, 0, -8 /)
      T(6,:) = (/ 1, 1, 1, -6, 1, 0, 0, -17 /)
      T(7,:) = (/ 1, 1, 1, -6, 0, 0, 0, -27 /)

      CALL THOMAS(7,3,3,T)
      WRITE(*,*) T(:,8)

      END PROGRAM EXAMPLE
! === T H O M A S   A L G O R I T H M   B A N D E D   M A T R I X ======
!     KL = Lower band: No. of sub-diagonals
!     KU = Upper band: No. of super-diagonals
!     If KL = KU = 1 then the solver works
!     similar to TDMA. The system of equations
!     has to be given to the solver in the
!     following compact form:
!     beginning from the left-most column
!     we fill T(:,j) with vectors containing
!     sub-diagonal, diagonal, super-diagonal
!     and finally the RHS (vector b) elements.
!     Example: N = 5, KL = 1, KU = 2
!     2  3  4  0  0 | 5
!     1  2  3  4  0 | 5
!     0  1  2  3  4 | 5
!     0  0  1  2  3 | 5
!     0  0  0  1  2 | 5
!     This system has to be rearranged to:
!     0  2  3  4 | 5
!     1  2  3  4 | 5
!     1  2  3  4 | 5
!     1  2  3  0 | 5
!     1  2  0  0 | 5

      SUBROUTINE THOMAS(N,KL,KU,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(N,KL+KU+2)

      DO K = 1, N-1
        NI = K + KL
        IF ( NI .GT. N ) NI = N
        DO I = K+1, NI
           U = T(I, K+KL-I+1) / T(K, KL+1)
           IF ( ABS(T(K, KL+1)) .LT. 1D-15 ) &
           WRITE(*,*) 'Check: diagonal element = 0'
           NJ = K + KL + KU - I + 1
           DO J = K+KL-I+2, NJ
              T(I,J) = T(I,J) - T(K, I+J-K) * U
           ENDDO
           T(I, KL+KU+2) = T(I, KL+KU+2) - T(K, KL+KU+2) * U
        ENDDO
      ENDDO

      DO I = N, 1, -1
         K = I + 1
         S = 0D0
         DO J = KL+2, KL+KU+1
            S = S + T(I,J) * T(K, KL+KU+2)
            K = K + 1
            IF ( K .GT. N ) EXIT
         ENDDO
         T(I, KL+KU+2) = ( T(I, KL+KU+2) - S ) / T(I, KL+1)
      ENDDO

      RETURN
      END SUBROUTINE
! ======================================================================
