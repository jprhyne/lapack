*> \brief \b DTRTRI_MOD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DTRTRI_MOD + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/DTRTRI_MOD.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/DTRTRI_MOD.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/DTRTRI_MOD.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRTRI_MOD( UPLO, DIAG, N, A, LDA, INFO )
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
*> DTRTRI_MOD computes the inverse of a real upper or lower triangular
*> matrix A.
*>
*> This is the Level 3 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  A is upper triangular;
*>          = 'L':  A is lower triangular.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          = 'N':  A is non-unit triangular;
*>          = 'U':  A is unit triangular.
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
*>          leading N-by-N upper triangular part of the array A contains
*>          the upper triangular matrix, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of the array A contains
*>          the lower triangular matrix, and the strictly upper
*>          triangular part of A is not referenced.  If DIAG = 'U', the
*>          diagonal elements of A are also not referenced and are
*>          assumed to be 1.
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
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*>               matrix is singular and its inverse can not be computed.
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
*> \ingroup trtri
*
*  =====================================================================
      SUBROUTINE DTRTRI_MOD( UPLO, DIAG, N, A, LDA, INFO )
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
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            K, J, JB, NB, NX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRI_MOD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine when to crossover to level2 version.
*
      !NX = ILAENV( 3, 'DTRTRI_MOD', UPLO // DIAG, N, -1, -1, -1 )
      NX = 1

      if( n.le.nx ) then
         call dtrti2_mod(uplo, diag, n, a, lda, info)
         return
      end if

      K = N/2

      ! X_{11} and X_{22} are computed with the same inputs everytime
      call DTRTRI_MOD(uplo, diag, n-k, a(k+1,k+1), lda, info)
      if( upper ) then
*
*        Break A apart as
*              |---------------|
*        A =   | A_{11} A_{12} |
*              | 0      A_{22} |
*              |---------------|
*
*              Then we want to overwrite A with the solution to
*
*              |---------------|    |---------------|   |-----|
*        A =   | A_{11} A_{12} | *  | X_{11} X_{12} | = | I 0 |
*              | 0      A_{22} |    | 0      X_{22} |   | 0 I |
*              |---------------|    |---------------|   |-----|
*
         ! X_{22} is computed above already

         ! if diag = 'u', then x is also unit
         call dtrmm('Right', 'Upper', 'No Transpose', diag, k, n-k,
     $         -one, a(k+1,k+1), lda, a(1,k+1), lda)

         ! finish computing X_{12}
         call dtrsm_mod( 'Left', 'Upper', 'No Transpose', diag, k, n-k,
     &         one, a, lda, a(1,k+1), lda )

         ! compute X_{11} afterwards
      else
*
*        Break A apart as
*              |---------------|
*        A =   | A_{11} 0      |
*              | A_{21} A_{22} |
*              |---------------|
*
*              Then we want to overwrite A with the solution to
*
*              |---------------|    |---------------|   |-----|
*        A =   | X_{11} 0      | *  | A_{11} 0      | = | I 0 |
*              | X_{21} X_{22} |    | A_{21} A_{22} |   | 0 I |
*              |---------------|    |---------------|   |-----|
*
         ! X_{22} is computed already
         ! if diag = 'u', then x is also unit
         call dtrmm('Left', 'Lower', 'No Transpose', diag, n-k, k,
     $      -one, a(k+1,k+1), lda, a(k+1,1), lda)

         ! finish computing X_{21}
         call dtrsmi_mod('Right', 'Lower', 'No Transpose', diag, n-k, k,
     $      one, a, lda, a(1,k+1), lda)
      end if

      call DTRTRI_MOD(uplo, diag, k, a, lda, info)

      END SUBROUTINE
