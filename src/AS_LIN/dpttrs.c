#include "f2c.h"

/* Subroutine */ int dpttrs_(integer *n, integer *nrhs, doublereal *d, 
	doublereal *e, doublereal *b, integer *ldb, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPTTRS solves a system of linear equations A * X = B with a   
    symmetric positive definite tridiagonal matrix A using the   
    factorization A = L*D*L**T or A = U**T*D*U computed by DPTTRF.   
    (The two forms are equivalent if A is real.)   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the tridiagonal matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the diagonal matrix D from the   
            factorization computed by DPTTRF.   

    E       (input) DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) off-diagonal elements of the unit bidiagonal factor 
  
            U or L from the factorization computed by DPTTRF.   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;
    /* Local variables */
    static integer i, j;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#define D(I) d[(I)-1]
#define E(I) e[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*ldb < max(1,*n)) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPTTRS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Solve A * X = B using the factorization A = L*D*L',   
       overwriting each right hand side vector with its solution. */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

/*        Solve L * x = b. */

	i__2 = *n;
	for (i = 2; i <= *n; ++i) {
	    B(i,j) -= B(i-1,j) * E(i - 1);
/* L10: */
	}

/*        Solve D * L' * x = b. */

	B(*n,j) /= D(*n);
	for (i = *n - 1; i >= 1; --i) {
	    B(i,j) = B(i,j) / D(i) - B(i+1,j) * E(i);
/* L20: */
	}
/* L30: */
    }

    return 0;

/*     End of DPTTRS */

} /* dpttrs_ */

