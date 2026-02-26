*> \brief \b ZTRTRSM solves the system XT = op(V) or TX = op(V) in place where all matrices are the same kind of triangular
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE ZTRTRSM(SIDE, UPLO, TRANS, DIAGT, DIAGV, N,
*    $            T, LDT, V, LDV)
*
*        .. Scalar Arguments ..
*        INTEGER           N, LDT, LDV
*        CHARACTER         SIDE, UPLO, TRANSV, DIAGT, DIAGV
*        ..
*        .. Array Arguments ..
*        COMPLEX*16        T(LDT,*), V(LDV,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRTRSM solves one of the following systems for X
*>
*>       T * X = op(V)
*>                      or
*>       X * T = op(V)
*> T and V are unit, or non-unit, upper or
*> lower triangular matrix, op(V) has the same shape as T , and op(V) is one of
*>
*>       op(V) = V      or       op(V) = V**T
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether T is on the left or right as
*>           follows:
*>
*>             SIDE = 'L' or 'l'    T*X = op(V)
*>
*>             SIDE = 'R' or 'r'    X*T = op(V)
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix T is an upper or
*>           lower triangular matrix as follows:
*>             UPLO = 'U' or 'u'    T is upper triangular
*>
*>             UPLO = 'L' or 'l'    T is lower triangular
*> \Endverbatim
*>
*> \param[in] TRANSV
*> \verbatim
*>          TRANSV is CHARACTER*1
*>           On entry, TRANSV specifies the form of op(V) to be used in
*>           the matrix multiplication as follows:
*>             TRANSV = 'N' or 'n'    op(V) = V
*>
*>             TRANSV = 'T' or 't'    op(V) = V**T
*>
*>             TRANSV = 'C' or 'c'    op(V) = V**T
*> \endverbatim
*>
*> \param[in] DIAGT
*> \verbatim
*>          DIAGT is CHARACTER*1
*>           On entry, DIAGT specifies whether or not T is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'      T is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'      T is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] DIAGV
*> \verbatim
*>          DIAGV is CHARACTER*1
*>           On entry, DIAGV specifies whether or not V is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'      V is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'      V is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of rows and columns of T.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in/out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension ( LDT, N )
*>           Before entry with UPLO = 'U' or 'u', the leading k-by-k
*>           upper triangular part of the array T must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           T is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l', the leading k-by-k
*>           lower triangular part of the array T must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           T is not referenced.
*>           Note that when  DIAGT = 'U' or 'u',  the diagonal elements of
*>           T  are not referenced either,  but are assumed to be  unity.
*>           On exit, T holds the array X
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>           On entry, LDT specifies the first dimension of T as declared
*>           in the calling (sub) program. LDT must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension ( LDV, N )
*>           Before entry with UPLO = 'U' or 'u', the leading k-by-k
*>           upper triangular part of the array op(V) must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           V is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l', the leading k-by-k
*>           lower triangular part of the array op(V) must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           V is not referenced.
*>           Note that when  DIAGV = 'U' or 'u',  the diagonal elements of
*>           V  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>           On entry, LDV specifies the first dimension of T as declared
*>           in the calling (sub) program. LDV must be at least max( 1, n ).
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
      SUBROUTINE ZTRTRSM(SIDE, UPLO, TRANS, DIAGT, DIAGV, N,
     $            T, LDT, V, LDV)
*
*        .. Scalar Arguments ..
         INTEGER           N, LDT, LDV
         CHARACTER         SIDE, UPLO, TRANS, DIAGT, DIAGV
*        ..
*        .. Array Arguments ..
         COMPLEX*16        T(LDT,*), V(LDV,*)
*        ..
*
*  =====================================================================
*
*        .. External Functions ..
         LOGICAL           LSAME
         EXTERNAL          LSAME
*        ..
*        .. Local Scalars ..
         INTEGER           I, J, K
         LOGICAL           TLEFT, TUPPER, VTRANS, VUNIT, TUNIT
         COMPLEX*16        TMP
*        ..
*        .. Local Parameters ..
         COMPLEX*16        ONE
         PARAMETER(ONE=(1.0D+0,0.0D+0))
*        ..
*
*        Beginning of Executable Statements
*
         TLEFT = LSAME(SIDE, 'L')
         TUPPER = LSAME(UPLO, 'U')
         VTRANS = LSAME(TRANS, 'T').OR.LSAME(TRANS, 'C')
         TUNIT = LSAME(DIAGT, 'U')
         VUNIT = LSAME(DIAGV, 'U')

         IF (TLEFT) THEN
            IF (TUPPER) THEN
               IF (VTRANS) THEN
                  IF (TUNIT) THEN
                     IF (VUNIT) THEN
*
*                       For each column in T
*
                        DO J = N, 1, -1
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = ONE
*
*                          compute all elements above
*
                           IF (J.GT.1) THEN
                              DO I = J-1, 1, -1
                                 TMP = V(J,I)
                                 DO K = J, I+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each column in T
*
                        DO J = N, 1, -1
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = V(J,J)
*
*                          Compute all elements above
*
                           IF (J.GT.1) THEN
                              DO I = J-1, 1, -1
                                 TMP = V(J,I)
                                 DO K = J, I+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
                     END IF
*                 .NOT.TUNIT
                  ELSE
                     IF (VUNIT) THEN
*
*                       For each column in T
*
                        DO J = N, 1, -1
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = ONE / T(J,J)
*
*                          Compute all elements above
*
                           IF (J.GT.1) THEN
                              DO I = J-1, 1, -1
                                 TMP = V(J,I)
                                 DO K = J, I+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(I,I)
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each column in T
*
                        DO J = N, 1, -1
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = V(J,J) / T(J,J)
*
*                          Compute all elements above
*
                           IF (J.GT.1) THEN
                              DO I = J-1, 1, -1
                                 TMP = V(J,I)
                                 DO K = J, I+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(I,I)
                              END DO
                           END IF
                        END DO
                     END IF
                  END IF
*              .NOT.VTRANS
               ELSE
                  IF (TUNIT) THEN
                     IF (VUNIT) THEN
*
*                       For each column in T
*
                        DO J = N, 1, -1
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = ONE
*
*                          Compute all elements above
*
                           IF (J.GT.1) THEN
                              DO I = J-1, 1, -1
                                 TMP = V(I,J)
                                 DO K = J, I+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each column in T
*
                        DO J = N, 1, -1
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = V(J,J)
*
*                          Compute all elements above
*
                           IF (J.GT.1) THEN
                              DO I = J-1, 1, -1
                                 TMP = V(I,J)
                                 DO K = J, I+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
                     END IF
*                 .NOT.TUNIT
                  ELSE
                     IF (VUNIT) THEN
*
*                       For each column in T
*
                        DO J = N, 1, -1
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = ONE / T(J,J)
*
*                          Compute all elements above
*
                           IF (J.GT.1) THEN
                              DO I = J-1, 1, -1
                                 TMP = V(I,J)
                                 DO K = J, I+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(I,I)
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each column in T
*
                        DO J = N, 1, -1
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = V(J,J) / T(J,J)
*
*                          Compute all elements above
*
                           IF (J.GT.1) THEN
                              DO I = J-1, 1, -1
                                 TMP = V(I,J)
                                 DO K = J, I+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(I,I)
                              END DO
                           END IF
                        END DO
                     END IF
                  END IF
               END IF
*           .NOT.TUPPER
            ELSE
               IF (VTRANS) THEN
                  IF (TUNIT) THEN
                     IF (VUNIT) THEN
*
*                       For each column in T
*
                        DO J = 1, N
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = ONE
*
*                          Compute all elements below
*
                           IF (J.LT.N) THEN
                              DO I = J+1, N
                                 TMP = V(J,I)
                                 DO K = J, I-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each column in T
*
                        DO J = 1, N
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = V(J,J)
*
*                          Compute all elements below
*
                           IF (J.LT.N) THEN
                              DO I = J+1, N
                                 TMP = V(J,I)
                                 DO K = J, I-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
                     END IF
*                 .NOT.TUNIT
                  ELSE
                     IF (VUNIT) THEN
*
*                       For each column in T
*
                        DO J = 1, N
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = ONE / T(J,J)
*
*                          Compute all elements below
*
                           IF (J.LT.N) THEN
                              DO I = J+1, N
                                 TMP = V(J,I)
                                 DO K = J, I-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(I,I)
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each column in T
*
                        DO J = 1, N
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = V(J,J) / T(J,J)
*
*                          Compute all elements below
*
                           IF (J.LT.N) THEN
                              DO I = J+1, N
                                 TMP = V(J,I)
                                 DO K = J, I-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(I,I)
                              END DO
                           END IF
                        END DO
                     END IF
                  END IF
*              .NOT.VTRANS
               ELSE
                  IF (TUNIT) THEN
                     IF (VUNIT) THEN
*
*                       For each column in T
*
                        DO J = 1, N
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = ONE
*
*                          Compute all elements below
*
                           IF (J.LT.N) THEN
                              DO I = J+1, N
                                 TMP = V(I,J)
                                 DO K = J, I-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each column in T
*
                        DO J = 1, N
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = V(J,J)
*
*                          Compute all elements below
*
                           IF (J.LT.N) THEN
                              DO I = J+1, N
                                 TMP = V(I,J)
                                 DO K = J, I-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
                     END IF
*                 .NOT.TUNIT
                  ELSE
                     IF (VUNIT) THEN
*
*                       For each column in T
*
                        DO J = 1, N
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = ONE / T(J,J)
*
*                          Compute all elements below
*
                           IF (J.LT.N) THEN
                              DO I = J+1, N
                                 TMP = V(I,J)
                                 DO K = J, I-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(I,I)
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each column in T
*
                        DO J = 1, N
*
*                          Compute from the bottom of the column to top
*                          Compute the element on the diagonal
*
                           T(J,J) = V(J,J) / T(J,J)
*
*                          Compute all elements below
*
                           IF (J.LT.N) THEN
                              DO I = J+1, N
                                 TMP = V(I,J)
                                 DO K = J, I-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(I,I)
                              END DO
                           END IF
                        END DO
                     END IF
                  END IF
               END IF
            END IF
*        .NOT.TLEFT
         ELSE
            IF (TUPPER) THEN
               IF (VTRANS) THEN
                  IF (TUNIT) THEN
                     IF (VUNIT) THEN
*
*                       For each row of T
*
                        DO I = 1, N
*
*                          Compute the leftmost component of the row
*
                           T(I,I) = ONE
*
*                          Compute all elements to the right (if any)
*
                           IF (I.LT.N) THEN
                              DO J = I+1, N
                                 TMP = V(J,I)
                                 DO K = I, J-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each row of T
*
                        DO I = 1, N
*
*                          Compute the leftmost component of the row
*
                           T(I,I) = V(I,I)
*
*                          Compute all elements to the right (if any)
*
                           IF (I.LT.N) THEN
                              DO J = I+1, N
                                 TMP = V(J,I)
                                 DO K = I, J-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
                     END IF
*                 .NOT.TUNIT
                  ELSE
                     IF (VUNIT) THEN
*
*                       For each row of T
*
                        DO I = 1, N
*
*                          Compute the leftmost component of the row
*
                           T(I,I) = ONE / T(I,I)
*
*                          Compute all elements to the right (if any)
*
                           IF (I.LT.N) THEN
                              DO J = I+1, N
                                 TMP = V(J,I)
                                 DO K = I, J-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(J,J)
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each row of T
*
                        DO I = 1, N
*
*                          Compute the leftmost component of the row
*
                           T(I,I) = V(I,I) / T(I,I)
*
*                          Compute all elements to the right (if any)
*
                           IF (I.LT.N) THEN
                              DO J = I+1, N
                                 TMP = V(J,I)
                                 DO K = I, J-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(J,J)
                              END DO
                           END IF
                        END DO
                     END IF
                  END IF
*              .NOT.VTRANS
               ELSE
                  IF (TUNIT) THEN
                     IF (VUNIT) THEN
*
*                       For each row of T
*
                        DO I = 1, N
*
*                          Compute the leftmost component of the row
*
                           T(I,I) = ONE
*
*                          Compute all elements to the right (if any)
*
                           IF (I.LT.N) THEN
                              DO J = I+1, N
                                 TMP = V(I,J)
                                 DO K = I, J-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each row of T
*
                        DO I = 1, N
*
*                          Compute the leftmost component of the row
*
                           T(I,I) = V(I,I)
*
*                          Compute all elements to the right (if any)
*
                           IF (I.LT.N) THEN
                              DO J = I+1, N
                                 TMP = V(I,J)
                                 DO K = I, J-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
                     END IF
*                 .NOT.TUNIT
                  ELSE
                     IF (VUNIT) THEN
*
*                       For each row of T
*
                        DO I = 1, N
*
*                          Compute the leftmost component of the row
*
                           T(I,I) = ONE / T(I,I)
*
*                          Compute all elements to the right (if any)
*
                           IF (I.LT.N) THEN
                              DO J = I+1, N
                                 TMP = V(I,J)
                                 DO K = I, J-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(J,J)
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each row of T
*
                        DO I = 1, N
*
*                          Compute the leftmost component of the row
*
                           T(I,I) = V(I,I) / T(I,I)
*
*                          Compute all elements to the right (if any)
*
                           IF (I.LT.N) THEN
                              DO J = I+1, N
                                 TMP = V(I,J)
                                 DO K = I, J-1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(J,J)
                              END DO
                           END IF
                        END DO
                     END IF
                  END IF
               END IF
*           .NOT.TUPPER
            ELSE
               IF (VTRANS) THEN
                  IF (TUNIT) THEN
                     IF (VUNIT) THEN
*
*                       For each row of T
*
                        DO I = N, 1, -1
*
*                          Compute the rightmost component of the row
*
                           T(I,I) = ONE
*
*                          If necessary, compute all elements to the left
*
                           IF (I.NE.1) THEN
                              DO J = I-1, 1, -1
                                 TMP = V(J,I)
                                 DO K = I, J+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each row of T
*
                        DO I = N, 1, -1
*
*                          Compute the rightmost component of the row
*
                           T(I,I) = V(I,I)
*
*                          If necessary, compute all elements to the left
*
                           IF (I.NE.1) THEN
                              DO J = I-1, 1, -1
                                 TMP = V(J,I)
                                 DO K = I, J+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
                     END IF
*                 .NOT.TUPPER
                  ELSE
                     IF (VUNIT) THEN
*
*                       For each row of T
*
                        DO I = N, 1, -1
*
*                          Compute the rightmost component of the row
*
                           T(I,I) = ONE / T(I,I)
*
*                          If necessary, compute all elements to the left
*
                           IF (I.NE.1) THEN
                              DO J = I-1, 1, -1
                                 TMP = V(J,I)
                                 DO K = I, J+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(J,J)
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each row of T
*
                        DO I = N, 1, -1
*
*                          Compute the rightmost component of the row
*
                           T(I,I) = V(I,I) / T(I,I)
*
*                          If necessary, compute all elements to the left
*
                           IF (I.NE.1) THEN
                              DO J = I-1, 1, -1
                                 TMP = V(J,I)
                                 DO K = I, J+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(J,J)
                              END DO
                           END IF
                        END DO
                     END IF
                  END IF
*              .NOT.VTRANS
               ELSE
                  IF (TUNIT) THEN
                     IF (VUNIT) THEN
*
*                       For each row of T
*
                        DO I = N, 1, -1
*
*                          Compute the rightmost component of the row
*
                           T(I,I) = ONE
*
*                          If necessary, compute all elements to the left
*
                           IF (I.NE.1) THEN
                              DO J = I-1, 1, -1
                                 TMP = V(I,J)
                                 DO K = I, J+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each row of T
*
                        DO I = N, 1, -1
*
*                          Compute the rightmost component of the row
*
                           T(I,I) = V(I,I)
*
*                          If necessary, compute all elements to the left
*
                           IF (I.NE.1) THEN
                              DO J = I-1, 1, -1
                                 TMP = V(I,J)
                                 DO K = I, J+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP
                              END DO
                           END IF
                        END DO
                     END IF
*                 .NOT.TUNIT
                  ELSE
                     IF (VUNIT) THEN
*
*                       For each row of T
*
                        DO I = N, 1, -1
*
*                          Compute the rightmost component of the row
*
                           T(I,I) = ONE / T(I,I)
*
*                          If necessary, compute all elements to the left
*
                           IF (I.NE.1) THEN
                              DO J = I-1, 1, -1
                                 TMP = V(I,J)
                                 DO K = I, J+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(J,J)
                              END DO
                           END IF
                        END DO
*                    .NOT.VUNIT
                     ELSE
*
*                       For each row of T
*
                        DO I = N, 1, -1
*
*                          Compute the rightmost component of the row
*
                           T(I,I) = V(I,I) / T(I,I)
*
*                          If necessary, compute all elements to the left
*
                           IF (I.NE.1) THEN
                              DO J = I-1, 1, -1
                                 TMP = V(I,J)
                                 DO K = I, J+1, -1
                                    TMP = TMP - T(I,K)*T(K,J)
                                 END DO
                                 T(I,J) = TMP / T(J,J)
                              END DO
                           END IF
                        END DO
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END SUBROUTINE
