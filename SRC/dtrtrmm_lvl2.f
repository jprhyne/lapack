*> \brief \b DTRTRMM_LVL2 computes an in place triangular-triangular matrix product
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE DTRTRMM_LVL2(SIDE, UPLO, TRANSV, DIAGT, DIAGV,
*    $                        N, ALPHA, T, LDT, V, LDV)
*
*        .. Scalar Arguments ..
*        INTEGER           N, LDT, LDV
*        CHARACTER         SIDE, UPLO, TRANSV, DIAGT, DIAGV
*        DOUBLE PRECISION  ALPHA
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
*> DTRTRMM_LVL2 performs one  of the matrix-matrix operations
*>
*>       T = \alpha op(V) * T
*>                      or
*>       T = \alpha T * op(V)
*> where \alpha is a scalar, T and V are unit, or non-unit, upper or
*> lower triangular matrix, and op(V) is one of
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
*>           On entry, SIDE specifies whether T multiplies op(V) from
*>           the left or right as follows:
*>
*>             SIDE = 'L' or 'l'    T = \alpha T * op(V)
*>
*>             SIDE = 'R' or 'r'    T = \alpha op(V) * T
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
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha. When alpha is
*>           zero then T and V are not referenced, and T and V need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in,out] T
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
*>           T  are not referenced either,  but are assumed to be  unit.
*>           But are set explicitly on exit
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
c     Cost: 2/21 * (3n^3 + 7n^2 - 10) + 2
      SUBROUTINE DTRTRMM_LVL2(SIDE, UPLO, TRANSV, DIAGT, DIAGV,
     $                        N, ALPHA, T, LDT, V, LDV)
*
*        .. Scalar Arguments ..
         INTEGER           N, LDT, LDV
         CHARACTER         SIDE, UPLO, TRANSV, DIAGT, DIAGV
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
*        .. External Subroutines ..
         EXTERNAL          DTRMM, DTRMMOOP, DLASET
*        ..
*        .. Local Scalars ..
         INTEGER           K, I
         CHARACTER         transv_2, uplov
         LOGICAL           TLEFT, TUPPER, VTRANS, VUNIT, TUNIT
*        ..
*        .. Local Parameters ..
         DOUBLE PRECISION ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0)
*        ..
*
*        Beginning of Executable Statements
*
*
*        Early Termination Criteria
*
         IF (ALPHA.EQ.ZERO) THEN
*
*           If ALPHA is 0, then we are just setting T to be the 0 matrix
*
            CALL DLASET(UPLO, N, N, ZERO, ZERO, T, LDT)
            RETURN
         END IF
*
*        Trivial Case
*
         IF(N.LE.0) THEN
            RETURN
         END IF
*
*        Convert our character flags into more understandible logical flags
*
         TUPPER = LSAME(UPLO,'U')
         TLEFT  = LSAME(SIDE,'L')
         VUNIT = LSAME(DIAGV, 'U')
         TUNIT = LSAME(DIAGT, 'U')
         VTRANS = LSAME(transv,'t').or.lsame(transv,'c')

         ! Must flip the transpose flag for V in some cases
         ! since dtrmv only accepts right multiplication from
         ! vectors.
         transv_2 = 't'
         if( vtrans ) then
            transv_2 = 'n'
         end if
         ! We need to know if v is upper or lower for passing into dtrmv
         ! if op(v) = v, then v is the same as T
         ! otherwise we must flip the value of uplo
         if( vtrans ) then
            if( tupper ) then
               uplov = 'l'
            else
               uplov = 'u'
            end if
         else
            uplov = uplo
         end if
         IF(TUPPER) THEN
*
*           T is upper triangular
*
            IF(TLEFT) THEN
*
*              Compute T = T*op(V)
*
               IF(TUNIT) THEN
                  ! DO later
               ELSE
                  do i = 1, n
                     ! scale the current row of t by alpha
                     call dscal( n-i+1, alpha, t(i,i), ldt )
                     call dtrmv( uplov, transv_2, diagv, n-i+1,
     $                  v(i,i), ldv, t(i,i), ldt )
                  end do
               END IF
            ELSE
*
*              Compute T = op(V)*T
*
               IF(VTRANS) THEN ! op(V) = V**T
               ELSE ! op(V) = V
               END IF
            END IF
         ELSE
*
*           T is lower triangular
*
            IF(TLEFT) THEN
*
*              Compute T = T*op(V)
*
               IF(VTRANS) THEN ! op(V) = V**T
               ELSE ! op(V) = V
               END IF
            ELSE
*
*              Compute T = op(V)*T
*
               IF(VTRANS) THEN ! op(V) = V**T
               ELSE ! op(V) = V
               END IF
            END IF
         END IF

      END SUBROUTINE
