*> \brief \b DTRTI2_MOD computes the inverse of a triangular matrix allowing for
*>    a zero diagonal (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DTRTI2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrti2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrti2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrti2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRTI2_mod computes the inverse of a real upper or lower triangular
*> matrix. It allows for a zero diagonal. If a diagonal is 0, then
*> that column is set to 0
*>
*> This is the Level 2 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the matrix A is upper or lower triangular.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          Specifies whether or not the matrix A is unit triangular.
*>          = 'N':  Non-unit triangular
*>          = 'U':  Unit triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the triangular matrix A.  If UPLO = 'U', the
*>          leading n by n upper triangular part of the array A contains
*>          the upper triangular matrix, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading n by n lower triangular part of the array A contains
*>          the lower triangular matrix, and the strictly upper
*>          triangular part of A is not referenced.  If DIAG = 'U', the
*>          diagonal elements of A are also not referenced and are
*>          assumed to be 1.
*>
*>          On exit, the (triangular) inverse of the original matrix, in
*>          the same storage format.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -k, the k-th argument had an illegal value
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
*> \ingroup trti2
*
*  =====================================================================
      SUBROUTINE DTRTI2_MOD( UPLO, DIAG, N, A, LDA, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UNIT, UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      UNIT = LSAME( DIAG, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UNIT .AND. .NOT.LSAME( DIAG, 'N' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO J = n, 1, -1
            ! Determine if we have a 0 diagonal.
            ! If not, we compute this column as normal
            ! Otherwise we set the current column to 0.
            ! Then all 0s should propogate to the right in the
            ! same column
            IF( UNIT ) THEN
               call dscal(j-1, -ONE, A(1,J), 1)
            ELSE
               ! This is the only case where we need to check
               ! for zero diagonals
               if( A(J,J).EQ.ZERO ) then
                  ! Set this column to 0
                  call dlaset('a', j, 1, zero, zero, a(1,j), lda)
               else
                  ! In this case we must manually invert the diagonal
                  ! element then propogate
                  A(J,J) = ONE / A(J,J)
                  call dscal(j-1, -A(J,J), A(1,J), 1)
               end if
            END IF
            ! Regardless of what we have done above, we want to replace
            ! the current column of A with the solution to
            ! A(1:j-1,1:j-1) x = A(1:j-1,j) and store this
            ! inside A(1:j-1,j) but potentially having some
            ! zeros on the diagonal of A, so we call our modified routine
            call DTRSM_MOD('Left', 'Upper', 'No Transpose', DIAG,
     $         j-1, 1, one, a, lda, a(1,j), lda)
         END DO
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO J = 1, n
            if( UNIT ) THEN
               call dscal(n-j, -ONE, A(J+1,J), 1)
            else
               if( A(j,j).eq.zero ) then
                  ! Set this column to 0
                  call dlaset('a', n-j, 1, zero, zero, a(j+1,j), lda)
                  !call dlaset('a', 1, j-1, zero, zero, a(j+1,1), lda)
               else
                  ! In this case we must manually invert the diagonal
                  ! element then propogate
                  A(J,J) = ONE / A(J,J)
                  call dscal(n-j, -A(J,J), A(j+1,J), 1)
               end if
            end if
            ! Regardless of what we have done above, we want to replace
            ! the current column of A with the solution to
            ! A(j+1:n,j+1:n) x = A(j+1:n,j) and store this
            ! inside A(j+1:n,j) but potentially having some
            ! zeros on the diagonal of A, so we call our modified routine
            call dtrsm_mod('Left', 'Lower', 'No Transpose', DIAG,
     $         n-j, 1, one, a(j+1,j+1), lda, a(j+1,j), lda)
         end do
      END IF
*
      RETURN
*
*     End of DTRTI2_MOD
*
      END
