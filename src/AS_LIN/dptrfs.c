#include "f2c.h"

/* Subroutine */ int dptrfs_(integer *n, integer *nrhs, doublereal *d, 
	doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer 
	*ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	 doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DPTRFS improves the computed solution to a system of linear   
    equations when the coefficient matrix is symmetric positive definite 
  
    and tridiagonal, and provides error bounds and backward error   
    estimates for the solution.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the tridiagonal matrix A.   

    E       (input) DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) subdiagonal elements of the tridiagonal matrix A.   

    DF      (input) DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the diagonal matrix D from the   
            factorization computed by DPTTRF.   

    EF      (input) DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) subdiagonal elements of the unit bidiagonal factor 
  
            L from the factorization computed by DPTTRF.   

    B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    X       (input/output) DOUBLE PRECISION array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by DPTTRS.   
            On exit, the improved solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= max(1,N).   

    FERR    (output) DOUBLE PRECISION array, dimension (NRHS)   
            The forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).   

    BERR    (output) DOUBLE PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Internal Parameters   
    ===================   

    ITMAX is the maximum number of steps of iterative refinement.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b11 = 1.;
    
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;
    /* Local variables */
    static doublereal safe1, safe2;
    static integer i, j;
    static doublereal s;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer count;
    static doublereal bi;
    extern doublereal dlamch_(char *);
    static doublereal cx, dx, ex;
    static integer ix;
    extern integer idamax_(integer *, doublereal *, integer *);
    static integer nz;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal lstres;
    extern /* Subroutine */ int dpttrs_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);
    static doublereal eps;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define DF(I) df[(I)-1]
#define EF(I) ef[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else if (*ldx < max(1,*n)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPTRFS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    FERR(j) = 0.;
	    BERR(j) = 0.;
/* L10: */
	}
	return 0;
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

    nz = 4;
    eps = dlamch_("Epsilon");
    safmin = dlamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;

/*     Do for each right hand side */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

	count = 1;
	lstres = 3.;
L20:

/*        Loop until stopping criterion is satisfied.   

          Compute residual R = B - A * X.  Also compute   
          abs(A)*abs(x) + abs(b) for use in the backward error bound. 
*/

	if (*n == 1) {
	    bi = B(1,j);
	    dx = D(1) * X(1,j);
	    WORK(*n + 1) = bi - dx;
	    WORK(1) = abs(bi) + abs(dx);
	} else {
	    bi = B(1,j);
	    dx = D(1) * X(1,j);
	    ex = E(1) * X(2,j);
	    WORK(*n + 1) = bi - dx - ex;
	    WORK(1) = abs(bi) + abs(dx) + abs(ex);
	    i__2 = *n - 1;
	    for (i = 2; i <= *n-1; ++i) {
		bi = B(i,j);
		cx = E(i - 1) * X(i-1,j);
		dx = D(i) * X(i,j);
		ex = E(i) * X(i+1,j);
		WORK(*n + i) = bi - cx - dx - ex;
		WORK(i) = abs(bi) + abs(cx) + abs(dx) + abs(ex);
/* L30: */
	    }
	    bi = B(*n,j);
	    cx = E(*n - 1) * X(*n-1,j);
	    dx = D(*n) * X(*n,j);
	    WORK(*n + *n) = bi - cx - dx;
	    WORK(*n) = abs(bi) + abs(cx) + abs(dx);
	}

/*        Compute componentwise relative backward error from formula 
  

          max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )   

          where abs(Z) is the componentwise absolute value of the matr
ix   
          or vector Z.  If the i-th component of the denominator is le
ss   
          than SAFE2, then SAFE1 is added to the i-th components of th
e   
          numerator and denominator before dividing. */

	s = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WORK(i) > safe2) {
/* Computing MAX */
		d__2 = s, d__3 = (d__1 = WORK(*n + i), abs(d__1)) / WORK(i);
		s = max(d__2,d__3);
	    } else {
/* Computing MAX */
		d__2 = s, d__3 = ((d__1 = WORK(*n + i), abs(d__1)) + safe1) / 
			(WORK(i) + safe1);
		s = max(d__2,d__3);
	    }
/* L40: */
	}
	BERR(j) = s;

/*        Test stopping criterion. Continue iterating if   
             1) The residual BERR(J) is larger than machine epsilon, a
nd   
             2) BERR(J) decreased by at least a factor of 2 during the
   
                last iteration, and   
             3) At most ITMAX iterations tried. */

	if (BERR(j) > eps && BERR(j) * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

	    dpttrs_(n, &c__1, &DF(1), &EF(1), &WORK(*n + 1), n, info);
	    daxpy_(n, &c_b11, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
		    ;
	    lstres = BERR(j);
	    ++count;
	    goto L20;
	}

/*        Bound error from formula   

          norm(X - XTRUE) / norm(X) .le. FERR =   
          norm( abs(inv(A))*   
             ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)   

          where   
            norm(Z) is the magnitude of the largest component of Z   
            inv(A) is the inverse of A   
            abs(Z) is the componentwise absolute value of the matrix o
r   
               vector Z   
            NZ is the maximum number of nonzeros in any row of A, plus
 1   
            EPS is machine epsilon   

          The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B)) 
  
          is incremented by SAFE1 if the i-th component of   
          abs(A)*abs(X) + abs(B) is less than SAFE2. */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WORK(i) > safe2) {
		WORK(i) = (d__1 = WORK(*n + i), abs(d__1)) + nz * eps * WORK(
			i);
	    } else {
		WORK(i) = (d__1 = WORK(*n + i), abs(d__1)) + nz * eps * WORK(
			i) + safe1;
	    }
/* L50: */
	}
	ix = idamax_(n, &WORK(1), &c__1);
	FERR(j) = WORK(ix);

/*        Estimate the norm of inv(A).   

          Solve M(A) * x = e, where M(A) = (m(i,j)) is given by   

             m(i,j) =  abs(A(i,j)), i = j,   
             m(i,j) = -abs(A(i,j)), i .ne. j,   

          and e = [ 1, 1, ..., 1 ]'.  Note M(A) = M(L)*D*M(L)'.   

          Solve M(L) * x = e. */

	WORK(1) = 1.;
	i__2 = *n;
	for (i = 2; i <= *n; ++i) {
	    WORK(i) = WORK(i - 1) * (d__1 = EF(i - 1), abs(d__1)) + 1.;
/* L60: */
	}

/*        Solve D * M(L)' * x = b. */

	WORK(*n) /= DF(*n);
	for (i = *n - 1; i >= 1; --i) {
	    WORK(i) = WORK(i) / DF(i) + WORK(i + 1) * (d__1 = EF(i), abs(d__1)
		    );
/* L70: */
	}

/*        Compute norm(inv(A)) = max(x(i)), 1<=i<=n. */

	ix = idamax_(n, &WORK(1), &c__1);
	FERR(j) *= (d__1 = WORK(ix), abs(d__1));

/*        Normalize error. */

	lstres = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    d__2 = lstres, d__3 = (d__1 = X(i,j), abs(d__1));
	    lstres = max(d__2,d__3);
/* L80: */
	}
	if (lstres != 0.) {
	    FERR(j) /= lstres;
	}

/* L90: */
    }

    return 0;

/*     End of DPTRFS */

} /* dptrfs_ */

