*> \brief \b DLUSM Solves a linear system where all matrices are triangular
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE DLUSM(SIDE, UPLOA, DIAGA, DIAGB, N, ALPHA,
*    $         A, LDA, WORK)
*
*     .. Scalar Arguments ..
*     INTEGER           N, LDA
*     CHARACTER         SIDE, UPLOA, DIAGA, DIAGB
*     DOUBLE PRECISION  ALPHA
*
*     .. Array Arguments ..
*     DOUBLE PRECISION  A(LDA,*), WORK(*)
*     ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLUSM Solves one of the following systems for X
*>
*>                A*X = \alpha*B
*>                      or
*>                X*A = \alpha*B
*>
*> where \alpha is a scalar, A, B are unit or non-unit triangular, exactly
*> one of which is Lower and Upper triangular, and at most one is non-unit
*>
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether A is on the left or right
*>           as follows:
*>
*>             SIDE = 'L' or 'l'    A*X = \alpha B
*>
*>             SIDE = 'R' or 'R'    X*A = \alpha B
*> \endverbatim
*>
*> \param[in] UPLOA
*> \verbatim
*>          UPLOA is CHARACTER*1
*>           On entry, UPLOA specifies whether A is upper or lower triangular
*>           as follows:
*>
*>              UPLOA = 'U' or 'u'      A is upper triangular
*>
*>              UPLOA = 'L' or 'l'      A is lower triangular
*> \endverbatim
*>
*> \param[in] DIAGA
*> \verbatim
*>          DIAGA is CHARACTER*1
*>           On entry, DIAGA specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAGL= 'U' or 'u'      A is assumed to be unit triangular.
*>
*>              DIAGL= 'N' or 'n'      A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] DIAGB
*> \verbatim
*>          DIAGB is CHARACTER*1
*>           On entry, DIAGB specifies whether or not B is unit triangular
*>           as follows:
*>
*>              DIAGU= 'U' or 'u'      B is assumed to be unit triangular.
*>
*>              DIAGU= 'N' or 'n'      B is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of rows and columns of A and B. N must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha. When alpha is
*>           zero then A is not referenced, and A need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, N ) where
*>           Before entry the leading n-by-n strictly upper triangular part of the array
*>           A must contain the upper triangular matrix U and the strictly lower triangular part of
*>           the leading n-by-n submatrix must contain the lower triangular matrix L.
*>           If DIAGL != 'U', then the diagonal is assumed to be part of L, and if
*>           DIAGU != 'U', then the diagonal is assumed to be part of U.
*>           Note: At most one of DIAGL and DIAGU can be not equal to 'U'.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*  =====================================================================
      SUBROUTINE DLUSM(SIDE, UPLOA, DIAGA, DIAGB, N, ALPHA,
     $         X, LDX, WORK)
*
*        .. Scalar Arguments ..
         INTEGER           N, LDX
         CHARACTER         SIDE, UPLOA, DIAGA, DIAGB
         DOUBLE PRECISION  ALPHA
*
*        .. Array Arguments ..
         DOUBLE PRECISION  X(LDX,*), WORK(*)
*        ..
*
*  =====================================================================
*
*        .. External Functions ..
         LOGICAL           LSAME
         DOUBLE PRECISION  DDOT
         EXTERNAL          LSAME, DDOT
*        ..
*        .. External Subroutines ..
         EXTERNAL          DGEMM, DTRMM, DLASET, XERBLA
*        ..
*        .. Local Scalars ..
         INTEGER           I, J
         LOGICAL           ALEFT, AUPPER, AUNIT, BUNIT
*        ..
*        .. Local Parameters ..
         DOUBLE PRECISION  ZERO
         PARAMETER(ZERO=0.0D+0)
*        ..
*        .. Beginning of executable statements
         ! Check our input flags
         ALEFT = LSAME(SIDE, 'L')
         AUPPER = LSAME(UPLOA, 'U')
         AUNIT = LSAME(DIAGA, 'U')
         BUNIT = LSAME(DIAGB, 'U')

         IF ((.NOT.AUNIT).AND.(.NOT.BUNIT)) THEN
            ! If we fail, then we report that the last diag argument is the
            ! error as we cannot know what the user meant
            CALL XERBLA('DLUSM', 4)
            RETURN
         END IF

         ! Case to ensure that we can assume n > 1 below
         IF (N.EQ.1) THEN
            ! in this case we are solving ax = \alpha b where all are scalars
            WORK(1) = ALPHA
            IF (.NOT.BUNIT) THEN
               WORK(1) = WORK(1) * X(1,1)
            END IF
            ! At this point, work = alpha*b
            ! divide by a if necessary
            IF (.NOT.AUNIT) THEN
               WORK(1) = WORK(1) / X(1,1)
            END IF
            X(1,1) = WORK(1)
            RETURN
         END IF

         IF (ALEFT) THEN
            ! This means we are solving A*X = \alpha * B
            ! We do this rowwise
            IF (AUPPER) THEN
               IF (BUNIT) THEN
                  IF (AUNIT) THEN
                     DO I = N, 1, -1
                        CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                        ! B is lower, so it starts at X(I,1) but ends at either
                        ! X(I,I) or X(I,I-1)
                        CALL DAXPY(I-1, ALPHA, X(I,1), LDX, WORK, 1)
                        WORK(I) = ALPHA
                        ! At this point, we have that work = alpha*B(I,:)
                        DO J = 1, N
                           ! subtract off A(I,I+1:N)*X(I+1:N,J)
                           WORK(J) = WORK(J) - DDOT(N-I, X(I,I+1), LDX,
     $                        X(I+1,J),1)
                        END DO
                        CALL DCOPY(N, WORK, 1, X(I,1), LDX)
                     END DO
                  ELSE ! not aunit
                     DO I = N, 1, -1
                        CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                        ! B is lower, so it starts at X(I,1) but ends at either
                        ! X(I,I) or X(I,I-1)
                        CALL DAXPY(I-1, ALPHA, X(I,1), LDX, WORK, 1)
                        WORK(I) = ALPHA
                        ! At this point, we have that work = alpha*B(I,:)
                        DO J = 1, N
                           ! subtract off A(I,I+1:N)*X(I+1:N,J)
                           WORK(J) = WORK(J) - DDOT(N-I, X(I,I+1), LDX,
     $                        X(I+1,J),1)
                        END DO
                        CALL DSCAL(N, 1/X(I,I), WORK, 1)
                        CALL DCOPY(N, WORK, 1, X(I,1), LDX)
                     END DO
                  END IF
               ELSE ! Not bunit, aunit must be true
                  DO I = N, 1, -1
                     CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                     ! B is lower, so it starts at X(I,1) but ends at either
                     ! X(I,I) or X(I,I-1)
                     CALL DAXPY(I, ALPHA, X(I,1), LDX, WORK, 1)
                     ! At this point, we have that work = alpha*B(I,:)
                     DO J = 1, N
                        ! subtract off A(I,I+1:N)*X(I+1:N,J)
                        WORK(J) = WORK(J) - DDOT(N-I, X(I,I+1), LDX,
     $                     X(I+1,J),1)
                     END DO
                     CALL DCOPY(N, WORK, 1, X(I,1), LDX)
                  END DO
               END IF
            ELSE
               IF (BUNIT) THEN
                  IF (AUNIT) THEN
                     DO I = 1, N
                        CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                        WORK(I) = ALPHA
                        CALL DAXPY(N-I, ALPHA, X(I,I+1), LDX,
     $                        WORK(I+1), 1)
                        ! At this point, we have that work = alpha*B(I,:)
                        DO J = 1, N
                           ! subtract off A(I,1:I-1)*X(1:I-1,J)
                           WORK(J) = WORK(J) - DDOT(I-1, X(I,1),
     $                        LDX, X(1,J),1)
                        END DO
                        CALL DCOPY(N, WORK, 1, X(I,1), LDX)
                     END DO
                  ELSE ! not aunit
                     DO I = 1, N
                        CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                        WORK(I) = ALPHA
                        CALL DAXPY(N-I, ALPHA, X(I,I+1), LDX,
     $                        WORK(I+1), 1)
                        ! At this point, we have that work = alpha*B(I,:)
                        DO J = 1, N
                           ! subtract off A(I,1:I-1)*X(1:I-1,J)
                           WORK(J) = WORK(J) - DDOT(I-1, X(I,1),
     $                        LDX, X(1,J), 1)
                        END DO
                        CALL DSCAL(N, 1/X(I,I), WORK, 1)
                        CALL DCOPY(N, WORK, 1, X(I,1), LDX)
                     END DO
                  END IF
               ELSE ! not bunit. We know that a must be unit
                  DO I = 1, N
                     CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                     CALL DAXPY(N-I+1, ALPHA, X(I,I), LDX,
     $                     WORK(I), 1)
                     ! At this point, we have that work = alpha*B(I,:)
                     DO J = 1, N
                        ! subtract off A(I,1:I-1)*X(1:I-1,J)
                        WORK(J) = WORK(J) - DDOT(I-1, X(I,1), LDX,
     $                     X(1,J),1)
                     END DO
                     CALL DCOPY(N, WORK, 1, X(I,1), LDX)
                  END DO
               END IF
            END IF
*        .NOT.ALEFT
         ELSE
            ! This means we are solving X*A = \alpha * B
            ! we do this columnwise
            IF (AUPPER) THEN
               ! A is upper triangular
               IF (BUNIT) THEN
                  IF (AUNIT) THEN
                     DO I = 1, N
                        CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                        WORK(I) = ALPHA
                        CALL DAXPY(N-I, ALPHA, X(I+1,I), 1,
     $                           WORK(I+1), 1)
                        ! now, work = alpha*b(:,i)
                        DO J = 1, N
                           ! subtract off x(j,1:i-1)*a(1:i-1,i)
                           WORK(J) = WORK(J) -
     $                        DDOT(I-1, X(J, 1), LDX, X(1,I), 1)
                        END DO
                        CALL DCOPY(N, WORK, 1, X(1,I), 1)
                     END DO
                  ELSE
                     DO I = 1, N
                        CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                        WORK(I) = ALPHA
                        CALL DAXPY(N-I, ALPHA, X(I+1,I), 1,
     $                        WORK(I+1), 1)
                        ! now, work = alpha*b(:,i)
                        DO J = 1, N
                           ! subtract off x(j,1:i-1)*a(1:i-1,i)
                           WORK(J) = WORK(J) -
     $                        DDOT(I-1, X(J,1), LDX, X(1,I), 1)
                        END DO
                        CALL DSCAL(N, 1/X(I,I), WORK, 1)
                        CALL DCOPY(N, WORK, 1, X(1,I), 1)
                     END DO
                  end if
               ELSE ! not bunit, so aunit must be true
                  DO I = 1, N
                     CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                     CALL DAXPY(N-I+1, ALPHA, X(I,I), 1, WORK(I), 1)
                     ! now, work = alpha*b(:,i)
                     DO J = 1, N
                        ! subtract off x(j,1:i-1)*a(1:i-1,i)
                        WORK(J) = WORK(J) - DDOT(I-1, X(J,1), LDX,
     $                     X(1,I), 1)
                     END DO
                     CALL DCOPY(N, WORK, 1, X(1,I), 1)
                  END DO
               END IF
            ELSE
               ! A is lower triangular
               IF (BUNIT) THEN
                  IF (AUNIT) THEN
                     DO I = N, 1, -1
                        CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                        CALL DAXPY(I-1, ALPHA, X(1,I), 1, WORK, 1)
                        WORK(I) = ALPHA
                        ! now, work = alpha*b(:,i)
                        DO J = 1, N
                           ! subtract off x(j,i+1:n)*a(i+1:n,i)
                           WORK(J) = WORK(J) -
     $                        DDOT(N-I, X(J,I+1), LDX, X(I+1,I), 1)
                        END DO
                        CALL DCOPY(N, WORK, 1, X(1,I), 1)
                     END DO
                  ELSE ! not aunit
                     DO I = N, 1, -1
                        CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                        CALL DAXPY(I-1, ALPHA, X(1,I), 1, WORK, 1)
                        WORK(I) = ALPHA
                        ! now, work = alpha*b(:,i)
                        DO J = 1, N
                           ! subtract off x(j,i+1:n)*a(i+1:n,i)
                           WORK(J) = WORK(J) -
     $                        DDOT(N-I, X(J,I+1), LDX, X(I+1,I), 1)
                        END DO
                        CALL DSCAL(N, 1/X(I,I), WORK, 1)
                        CALL DCOPY(N, WORK, 1, X(1,I), 1)
                     END DO
                  END IF
               ELSE ! not bunit, so aunit must be true
                  DO I = N, 1, -1
                     CALL DLASET('A', N, 1, ZERO, ZERO, WORK, 1)
                     CALL DAXPY(I, ALPHA, X(1,I), 1, WORK, 1)
                     ! now, work = alpha*b(:,i)
                     DO J = 1, N
                        ! subtract off x(j,i+1:n)*a(i+1:n,i)
                        WORK(J) = WORK(J) -
     $                     DDOT(N-I, X(J,I+1), LDX, X(I+1,I), 1)
                     END DO
                     CALL DCOPY(N, WORK, 1, X(1,I), 1)
                  END DO
               END IF
            END IF
         END IF
      END SUBROUTINE
