*> \brief \b DLALL2 computes the product LLH or UHU, where U and L are upper or lower triangular matrices (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download dlall2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlall2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlall2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlall2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLALL2( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
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
*> DLALL2 computes the product L * L**H or U**H * U, where the triangular
*> factor U or L is stored in the upper or lower triangular part of
*> the array A.
*>
*> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
*> overwriting the factor U in A.
*> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
*> overwriting the factor L in A.
*>
*> This is the unblocked form of the algorithm, calling Level 2 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the triangular factor stored in the array A
*>          is upper or lower triangular:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the triangular factor U or L.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the triangular factor U or L.
*>          On exit, if UPLO = 'U', the upper triangle of A is
*>          overwritten with the upper triangle of the product U**h * U;
*>          if UPLO = 'L', the lower triangle of A is overwritten with
*>          the lower triangle of the product L * L**H.
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
*> \ingroup lall2
*
*  =====================================================================
      SUBROUTINE DLALL2( UPLO, N, A, LDA, INFO )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I,J
      DOUBLE PRECISION   AII,AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      !DOUBLE PRECISION   DDOT
      !EXTERNAL           LSAME, DDOT
*     ..
*     .. External Subroutines ..
      !EXTERNAL           DGEMV, DSCAL, XERBLA
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
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLALL2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( UPPER ) THEN
         ! Compute and store ut(U**H U)
         DO J = N, 1, -1
            ! Complex will want to store the conjugate
            AJJ = A(J,J)
            IF( J.EQ.1 ) THEN
               CALL DSCAL(N, AJJ, A, LDA)
            ELSE
               ! complex version will need {c,z}gemm since the vector x is
               ! conjugated (or a {c,z}gemcv-like statement)
               CALL DGEMV('Transpose', J-1, N-J+1, ONE,
     $            A(1,J), LDA, A(1,J), 1, AJJ, A(J,J), LDA)
            END IF
         END DO
      ELSE
         ! compute and store lt(LL**H)
         DO J = N, 1, -1
            ! complex will want to store the conjugate
            AJJ = A(J,J)
            IF( J.EQ.1 ) THEN
               ! n is simplified from n-j+1 for readability
               ! a=a(1,1) is simplified from a(j,j) also for readability
               CALL DSCAL(N, AJJ, A, 1)
            ELSE
               ! complex version will need {c,z}gemm since the vector x is
               ! conjugated (or a {c,z}gemcv-like statement)
               CALL DGEMV('No Transpose', N-J+1, J-1, ONE,
     $            A(J,1), LDA, A(J,1), LDA, AJJ, A(J,J), 1)
            END IF
         END DO
      END IF
      END SUBROUTINE
