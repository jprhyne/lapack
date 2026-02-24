*> \brief \b DTRTRMS solves the system XT = op(V) or TX = op(V) in place where all matrices are the same kind of triangular
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     RECURSIVE SUBROUTINE DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV,
*    $            N, ALPHA, T, LDT, V, LDV)
*
*        .. Scalar Arguments ..
*        INTEGER           N, LDT, LDV
*        CHARACTER         SIDE, UPLO, TRANSV, DIAGT, DIAGV
*        ..
*        .. Array Arguments ..
*        DOUBLE PRECISION  T(LDT,*), V(LDV,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRTRMS solves one of the following systems for X
*>
*>       T * X = alpha * op(V)
*>                      or
*>       X * T = alpha * op(V)
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
*> \endverbatim
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
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION
*>             On entry, ALPHA specifies the scalar alpha
*> \endverbatim
*>
*> \param[in/out] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension ( LDT, N )
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
*>          V is DOUBLE PRECISION array, dimension ( LDV, N )
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
      RECURSIVE SUBROUTINE DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV,
     $            N, ALPHA, T, LDT, V, LDV)
*
*        .. Scalar Arguments ..
         INTEGER           N, LDT, LDV
         CHARACTER         SIDE, UPLO, TRANS, DIAGT, DIAGV
         DOUBLE PRECISION  ALPHA
*        ..
*        .. Array Arguments ..
         DOUBLE PRECISION  T(LDT,*), V(LDV,*)
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
*        ..
*        .. Local Parameters ..
         DOUBLE PRECISION  NEG_ONE, ONE
         PARAMETER(NEG_ONE=-1.0D+0, ONE=1.0D+0)
*        ..
*
*        Beginning of Executable Statements
*
         TLEFT = LSAME(SIDE, 'L')
         TUPPER = LSAME(UPLO, 'U')
         VTRANS = LSAME(TRANS, 'T').OR.LSAME(TRANS, 'C')
         TUNIT = LSAME(DIAGT, 'U')
         VUNIT = LSAME(DIAGV, 'U')
*
*        Terminating case
*
         IF (N.EQ.1) THEN
            CALL DTRTRMS_LVL2(SIDE, UPLO, TRANS, DIAGT, DIAGV,
     $            N, ALPHA, T, LDT, V, LDV)
            RETURN
         ELSE IF(N.LE.0) THEN
            RETURN
         END IF
*
*        Recursive case
*
         K = N/2
         IF (TUPPER) THEN
*
*           T is upper triangular
*
            IF (TLEFT) THEN
*
*              Solve T*X = alpha * op(V) for X overwritting T
*
               IF (VTRANS) THEN
*
*                 We are computing X such that T*X = V**T, which we break down as follows
*                 |--------------|   |--------------|            |--------------------|
*                 |T_{11}  T_{12}|   |X_{11}  X_{12}|            |V_{11}**T  V_{21}**T|
*                 |0       T_{22}| * |0       X_{22}|  = alpha * |0          V_{22}**T|
*                 |--------------|   |--------------|            |--------------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}    T_{12}\in\R^{k\times n-k}
*                                            T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}
*                 V_{21}\in\R^{n-k\times k}  V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get the following systems we need to solve
*
*                 T_{11}*X_{11}                 = alpha V_{11}^\top
*                 T_{11}*X_{12} + T_{12}*X_{22} = alpha V_{21}^\top
*                 T_{22}*X_{22}                 = alpha V_{22}^\top
*
*                 Computing X_{11} and X_{22} are just recursive calls to this
*                 routine, Which we must do in reverse order since we cannot
*                 overwrite T_{11} until after we compute X_{12}, which we need
*                 X_{22} first.
*
*                 Solve T_{22}*X_{22} = alpha V_{22}^\top
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, N-K,
     $                  ALPHA, T(K+1,K+1), LDT, V(K+1,K+1), LDV)
*
*                 We solve the system
*
*                 T_{11}*X_{12} + T_{12}*X_{22} = alpha V_{21}^\top
*
*                 Our final result will be stored in T_{12}, so we can use
*                 this as a temporary workspace
*
*                 T_{12} = -T_{12}*X_{22}        (TRMM)
*
                  CALL DTRMM('Right', UPLO, 'No Transpose', 'Non-unit',
     $                  K, N-K, NEG_ONE, T(K+1,K+1), LDT, T(1,K+1), LDT)
*
*                 T_{12} = T_{12} + alpha V_{21}^\top
*
                  DO I = 1, K
                     DO J = K+1, N
                        T(I,J) = T(I,J) + ALPHA * V(J,I)
                     END DO
                  END DO
*
*                 Solve T_{11}*X_{12} = T_{12}  (TRSM)
*
                  CALL DTRSM('Left', UPLO, 'No Transpose', DIAGT,
     $                  K, N-K, ONE, T, LDT, T(1,K+1), LDT)
*
*                 Solve T_{11}*X_{11} = alpha V_{11}^\top  (TRTRMS)
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, K,
     $                  ALPHA, T, LDT, V, LDV)
*              .NOT.VTRANS
               ELSE
*
*                 We are computing X such that T*X = alpha V, which we break down as follows
*                 |--------------|   |--------------|            |--------------|
*                 |T_{11}  T_{12}|   |X_{11}  X_{12}|            |V_{11}  V_{12}|
*                 |0       T_{22}| * |0       X_{22}|  = alpha * |0       V_{22}|
*                 |--------------|   |--------------|            |--------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}    T_{12}\in\R^{k\times n-k}
*                                            T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}    V_{12}\in\R^{k\times n-k}
*                                            V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get the following systems we need to solve
*
*                 T_{11}*X_{11}                 = alpha V_{11}
*                 T_{11}*X_{12} + T_{12}*X_{22} = alpha V_{12}
*                 T_{22}*X_{22}                 = alpha V_{22}
*
*                 Computing X_{11} and X_{22} are just recursive calls to this
*                 routine, Which we must do in reverse order since we cannot
*                 overwrite T_{11} until after we compute X_{12}, which we need
*                 X_{22} first.
*
*                 Solve T_{22}*X_{22} = alpha V_{22}
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, N-K,
     $                  ALPHA, T(K+1,K+1), LDT, V(K+1,K+1), LDV)
*
*                 We solve the system
*
*                 T_{11}*X_{12} + T_{12}*X_{22} = alpha V_{12}
*
*                 Our final result will be stored in T_{12}, so we can use
*                 this as a temporary workspace
*
*                 T_{12} = -T_{12}*X_{22}        (TRMM)
*
                  CALL DTRMM('Right', UPLO, 'No Transpose', 'Non-unit',
     $                  K, N-K, NEG_ONE, T(K+1,K+1), LDT, T(1,K+1), LDT)
*
*                 T_{12} = T_{12} + alpha V_{12}
*
                  DO I = 1, K
                     DO J = K+1, N
                        T(I,J) = T(I,J) + ALPHA * V(I,J)
                     END DO
                  END DO
*
*                 Solve T_{11}*X_{12} = T_{12}  (TRSM)
*
                  CALL DTRSM('Left', UPLO, 'No Transpose', DIAGT,
     $                  K, N-K, ONE, T, LDT, T(1,K+1), LDT)
*
*                 Solve T_{11}*X_{11} = alpha V_{11}  (TRTRMS)
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, K,
     $                  ALPHA, T, LDT, V, LDV)
               END IF
            ELSE
*
*              Solve X*T = alpha op(V) for X overwriting T
*
               IF (VTRANS) THEN
*
*                 We are computing X such that X*T = alpha V**T
*                 |--------------|   |--------------|           |--------------------|
*                 |X_{11}  X_{12}|   |T_{11}  T_{12}|           |V_{11}**T  V_{21}**T|
*                 |0       X_{22}| * |0       T_{22}| = alpha * |0          V_{22}**T|
*                 |--------------|   |--------------|           |--------------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}    T_{12}\in\R^{k\times n-k}
*                                            T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}
*                 V_{21}\in\R^{n-k\times k}  V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we need to solve the following systems
*
*                 X_{11}*T_{11}                 = alpha V_{11}**T
*                 X_{11}*T_{12} + X_{12}*T_{22} = alpha V_{21}**T
*                 X_{22}*T_{22}                 = alpha V_{22}**T
*
*                 we compute X_{11} since this will be used to help compute
*                 X_{12}.
*
*                 Solve X_{11}*T_{11} = alpha V_{11}**T
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, K,
     $                  ALPHA, T, LDT, V, LDV)
*
*                 We solve the system
*
*                 X_{11}*T_{12} + X_{12}*T_{22} = alpha V_{21}**T
*
*                 Our final result will be stored in T_{12}, so we can use
*                 this as a temporary workspace
*
*                 T_{12} = -X_{11}*T_{12}        (TRMM)
*
                  CALL DTRMM('Left', UPLO, 'No Transpose',
     $                  'Non-unit', K, N-K, NEG_ONE, T, LDT, T(1,K+1),
     $                  LDT)
*
*                 T_{12} = T_{12} + alpha V_{21}^\top
*
                  DO I = 1, K
                     DO J = K+1, N
                        T(I,J) = T(I,J) + ALPHA * V(J,I)
                     END DO
                  END DO
*
*                 Solve X_{12}*T_{22} = T_{12}  (TRSM)
*
                  CALL DTRSM('Right', UPLO, 'No Transpose', DIAGT,
     $                  K, N-K, ONE, T(K+1,K+1), LDT, T(1,K+1), LDT)
*
*                 Solve T_{22}*X_{22} = alpha V_{22}^\top
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, N-K,
     $                  ALPHA, T(K+1,K+1), LDT, V(K+1,K+1), LDV)
               ELSE
*
*                 We are computing X such that X*T = alpha V
*                 |--------------|   |--------------|           |--------------|
*                 |X_{11}  X_{12}|   |T_{11}  T_{12}|           |V_{11}  V_{21}|
*                 |0       X_{22}| * |0       T_{22}| = alpha * |0       V_{22}|
*                 |--------------|   |--------------|           |--------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}    T_{12}\in\R^{k\times n-k}
*                                            T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}    V_{12}\in\R^{k\times n-k}
*                                            V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we need to solve the following systems
*
*                 X_{11}*T_{11}                 = alpha V_{11}
*                 X_{11}*T_{12} + X_{12}*T_{22} = alpha V_{12}
*                 X_{22}*T_{22}                 = alpha V_{22}
*
*                 we compute X_{11} since this will be used to help compute
*                 X_{12}.
*
*                 Solve X_{11}*T_{11} = V_{11}
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, K,
     $                  ALPHA, T, LDT, V, LDV)
*
*                 We solve the system
*
*                 X_{11}*T_{12} + X_{12}*T_{22} = alpha V_{12}
*
*                 Our final result will be stored in T_{12}, so we can use
*                 this as a temporary workspace
*
*                 T_{12} = -X_{11}*T_{12}        (TRMM)
*
                  CALL DTRMM('Left', UPLO, 'No Transpose',
     $                  'Non-unit', K, N-K, NEG_ONE, T, LDT, T(1,K+1),
     $                  LDT)
*
*                 T_{12} = T_{12} + alpha V_{12}
*
                  DO I = 1, K
                     DO J = K+1, N
                        T(I,J) = T(I,J) + ALPHA * V(I,J)
                     END DO
                  END DO
*
*                 Solve X_{12}*T_{22} = T_{12}  (TRSM)
*
                  CALL DTRSM('Right', UPLO, 'No Transpose', DIAGT,
     $                  K, N-K, ONE, T(K+1,K+1), LDT, T(1,K+1), LDT)
*
*                 Solve T_{22}*X_{22} = alpha V_{22}^\top
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, N-K,
     $                  ALPHA, T(K+1,K+1), LDT, V(K+1,K+1), LDV)
               END IF
            END IF
*        .NOT.TUPPER
         ELSE
*
*           T is lower triangular
*
            IF (TLEFT) THEN
*
*              Solve T*X = alpha op(V) for X overwritting T
*
               IF (VTRANS) THEN
*
*                 We are computing X such that T*X = V**T, which we break down as follows
*                 |--------------|   |--------------|            |--------------------|
*                 |T_{11}  0     |   |X_{11}  0     |            |V_{11}**T  0        |
*                 |T_{21}  T_{22}| * |X_{21}  X_{22}|  = alpha * |V_{12}**T  V_{22}**T|
*                 |--------------|   |--------------|            |--------------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}    V_{12}\in\R^{k\times n-k}
*                                            V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get the following systems we need to solve
*
*                 T_{11}*X_{11}                 = alpha V_{11}**T
*                 T_{21}*X_{11} + T_{22}*X_{21} = alpha V_{12}**T
*                 T_{22}*X_{22}                 = alpha V_{22}**T
*
*                 we compute X_{11} since this will be used to help compute
*                 X_{21}.
*
*                 Solve T_{11}*X_{11} = alpha V_{11}**T
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, K,
     $                  ALPHA, T, LDT, V, LDV)
*
*                 We solve the system
*
*                 T_{21}*X_{11} + T_{22}*X_{21} = alpha V_{12}**T
*
*                 Our final result will be stored in T_{21}, so we can use
*                 this as a temporary workspace
*
*                 T_{21} = -T_{21}*X_{11} (TRMM)
*
                  CALL DTRMM('Right', UPLO, 'No Transpose',
     $                  'Non-unit', N-K, K, NEG_ONE, T, LDT, T(K+1,1),
     $                  LDT)
*
*                 T_{21} = T_{21} + alpha V_{12}**T
*
                  DO I = K+1, N
                     DO J = 1, K
                        T(I,J) = T(I,J) + ALPHA * V(J,I)
                     END DO
                  END DO
*
*                 Solve T_{22}*X_{21} = T_{21}  (TRSM)
*
                  CALL DTRSM('Left', UPLO, 'No Transpose', DIAGT,
     $                  N-K, K, ONE, T(K+1,K+1), LDT, T(K+1,1), LDT)
*
*                 Solve T_{22}*X_{22} = alpha V_{22}**T
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, N-K,
     $                  ALPHA, T(K+1,K+1), LDT, V(K+1,K+1), LDV)
*              .NOT.VTRANS
               ELSE
*
*                 We are computing X such that T*X = alpha V**T, which we break down as follows
*                 |--------------|   |--------------|            |--------------|
*                 |T_{11}  0     |   |X_{11}  0     |            |V_{11}  0     |
*                 |T_{21}  T_{22}| * |X_{21}  X_{22}|  = alpha * |V_{21}  V_{22}|
*                 |--------------|   |--------------|            |--------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get the following systems we need to solve
*
*                 T_{11}*X_{11}                 = alpha V_{11}
*                 T_{21}*X_{11} + T_{22}*X_{21} = alpha V_{21}
*                 T_{22}*X_{22}                 = alpha V_{22}
*
*                 we compute X_{11} since this will be used to help compute
*                 X_{21}.
*
*                 Solve T_{11}*X_{11} = alpha V_{11}
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, K,
     $                  ALPHA, T, LDT, V, LDV)
*
*                 We solve the system
*
*                 T_{21}*X_{11} + T_{22}*X_{21} = alpha V_{21}
*
*                 Our final result will be stored in T_{21}, so we can use
*                 this as a temporary workspace
*
*                 T_{21} = -T_{21}*X_{11} (TRMM)
*
                  CALL DTRMM('Right', UPLO, 'No Transpose',
     $                  'Non-unit', N-K, K, NEG_ONE, T, LDT, T(K+1,1),
     $                  LDT)
*
*                 T_{21} = T_{21} + alpha V_{21}
*
                  DO I = K+1, N
                     DO J = 1, K
                        T(I,J) = T(I,J) + ALPHA * V(I,J)
                     END DO
                  END DO
*
*                 Solve T_{22}*X_{21} = T_{21}  (TRSM)
*
                  CALL DTRSM('Left', UPLO, 'No Transpose', DIAGT,
     $                  N-K, K, ONE, T(K+1,K+1), LDT, T(K+1,1), LDT)
*
*                 Solve T_{22}*X_{22} = alpha V_{22}
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, N-K,
     $                  ALPHA, T(K+1,K+1), LDT, V(K+1,K+1), LDV)
               END IF
            ELSE
*
*              Compute X such that X*T = alpha op(V)
*
               IF (VTRANS) THEN
*
*                 We are computing X such that X*T = alpha V**T, which we break down as follows
*                 |--------------|   |--------------|            |-------------------|
*                 |X_{11}  0     |   |T_{11}  0     |            |V_{11}**T 0        |
*                 |X_{21}  X_{22}| * |T_{21}  T_{22}|  = alpha * |V_{12}**T V_{22}**T|
*                 |--------------|   |--------------|            |-------------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get the following systems we need to solve
*
*                 X_{11}*T_{11}                 = alpha V_{11}**T
*                 X_{21}*T_{11} + X_{22}*T_{21} = alpha V_{12}**T
*                 X_{22}*T_{22}                 = alpha V_{22}**T
*
*                 We compute X_{22} since this will be used to help compute
*                 X_{21}.
*
*                 Solve X_{22}*T_{22} = alpha V_{22}**T
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, N-K,
     $                  ALPHA, T(K+1,K+1), LDT, V(K+1,K+1), LDV)
*
*                 We solve the system
*
*                 X_{21}*T_{11} + X_{22}*T_{21} = alpha V_{12}**T
*
*                 Our final result will be stored in T_{21}, so we can use
*                 this as a temporary workspace
*
*                 T_{21} = -X_{22}*T_{21} (TRMM)
*
                  CALL DTRMM('Left', UPLO, 'No Transpose',
     $                  'Non-Unit', N-K, K, NEG_ONE, T(K+1,K+1), LDT,
     $                  T(K+1,1), LDT)
*
*                 T_{21} = T_{21} + alpha V_{12}**T
*
                  DO I = K+1, N
                     DO J = 1, K
                        T(I,J) = T(I,J) + ALPHA * V(J,I)
                     END DO
                  END DO
*
*                 Solve X_{21}*T_{11} = T_{21}  (TRSM)
*
                  CALL DTRSM('Right', UPLO, 'No Transpose', DIAGT,
     $                  N-K, K, ONE, T, LDT, T(K+1,1), LDT)
*
*                 Solve X_{11}*T_{11} = alpha V_{11}**T
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, K,
     $                  ALPHA, T, LDT, V, LDV)
               ELSE
*
*                 We are computing X such that X*T = alpha V**T, which we break down as follows
*                 |--------------|   |--------------|            |-------------|
*                 |X_{11}  0     |   |T_{11}  0     |            |V_{11} 0     |
*                 |X_{21}  X_{22}| * |T_{21}  T_{22}|  = alpha * |V_{21} V_{22}|
*                 |--------------|   |--------------|            |-------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get the following systems we need to solve
*
*                 X_{11}*T_{11}                 = alpha V_{11}
*                 X_{21}*T_{11} + X_{22}*T_{21} = alpha V_{21}
*                 X_{22}*T_{22}                 = alpha V_{22}
*
*                 We compute X_{22} since this will be used to help compute
*                 X_{21}.
*
*                 Solve X_{22}*T_{22} = alpha V_{22}
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, N-K,
     $                  ALPHA, T(K+1,K+1), LDT, V(K+1,K+1), LDV)
*
*                 We solve the system
*
*                 X_{21}*T_{11} + X_{22}*T_{21} = alpha V_{21}
*
*                 Our final result will be stored in T_{21}, so we can use
*                 this as a temporary workspace
*
*                 T_{21} = X_{22}*T_{21} (TRMM)
*
                  CALL DTRMM('Left', UPLO, 'No Transpose',
     $                  'Non-Unit', N-K, K, NEG_ONE, T(K+1,K+1), LDT,
     $                  T(K+1,1), LDT)
*
*                 T_{21} = T_{21} + alpha V_{21}
*
                  DO I = K+1, N
                     DO J = 1, K
                        T(I,J) = T(I,J) + ALPHA * V(I,J)
                     END DO
                  END DO
*
*                 Solve X_{21}*T_{11} = T_{21}  (TRSM)
*
                  CALL DTRSM('Right', UPLO, 'No Transpose', DIAGT,
     $                  N-K, K, ONE, T, LDT, T(K+1,1), LDT)
*
*                 Solve X_{11}*T_{11} = alpha V_{11}
*
                  CALL DTRTRMS(SIDE, UPLO, TRANS, DIAGT, DIAGV, K,
     $                  ALPHA, T, LDT, V, LDV)
               END IF
            END IF
         END IF
      END SUBROUTINE
