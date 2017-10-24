/*==========================================================
 * adi_tri.c
 *
 * Solves the equation TX + XT' = F, where T is a special
 * tridiagonal matrix.
 *
 * Inputs:
 *     d = main diagonal of T
 *     u = upper diagonal of T
 *     l = lower diagonal of T
 *     p = shift parameters
 *     F = right-hand side matrix
 *
 * Outputs:
 *     X, where X is the result of
 *
 *     X = zeros(n,n);
 *     for j = 1:numel(p)
 *         X = (T+p(j)*I) \ ( F' - ((T-p(j)*I)*X)' );
 *         X = (T+p(j)*I) \ ( F  - ((T-p(j)*I)*X)' );
 *     end
 *
 * The calling syntax is:
 *
 *     X = adi_tri( d, u, l, p, F )
 *
 *========================================================*/

#include "mex.h"

void adi_tri(double *d, double *u, double *l, mwSize n, mwSize nshifts, double *p, double *F, double *X)
{
    int i, j, k;
    double s, *w, *v, *Y;

    w = (double *)mxMalloc(n * sizeof(double));
    v = (double *)mxMalloc(n * sizeof(double));
    Y = (double *)mxMalloc(n * n * sizeof(double));
    
    /* Optimization: Assuming X = 0 initially lets us simplify the first
     *               iteration without explicitly setting X to zero.
     *
     * See main loop for notes on the optimizations performed here.
     */

    s = p[0];

    /* Solve:
     *  
     *     X = (T+p(j)*I) \ F'
     */

    w[0] = d[0] + s;
    for (j = 0; j < n; ++j) {
        Y[j] = F[j] / w[0];
    }

    for (i = 1; i < n; ++i) {
        v[i-1] = u[i-1] / w[i-1];
        w[i] = (d[i]+s) - l[i-1]*v[i-1];
        Y[i*n] = ((F[i*n] - (d[0]-s)*X[i*n] - u[0]*X[i*n+1]) - l[i-1]*Y[(i-1)*n]) / w[i];
        for (j = 1; j < n-1; ++j) {
            Y[i*n+j] = ((F[i*n+j] - l[j-1]*X[i*n+(j-1)] - (d[j]-s)*X[i*n+j] - u[j]*X[i*n+(j+1)])
                        - l[i-1]*Y[(i-1)*n+j]) / w[i];
        }
        Y[i*n+(n-1)] = ((F[i*n+(n-1)] - l[n-2]*X[i*n+(n-2)] - (d[n-1]-s)*X[i*n+(n-1)]) - l[i-1]*Y[(i-1)*n+(n-1)]) / w[i];
    }

    /* Back substitution */

    for (i = n-2; i >= 0; --i) {
        for (j = 0; j < n; ++j) {
            Y[i*n+j] = Y[i*n+j] - v[i]*Y[(i+1)*n+j];
        }
    }

    /* Solve:
     * 
     *    X = (T+p(j)*I) \ ( F - ((T-p(j)*I)*X)' )
     */

    /* Forward pass */

    i = 0;
    X[i*n] = (F[i*n] - (d[i]-s)*Y[i*n] - u[i]*Y[(i+1)*n]) / w[0];
    for (j = 1; j < n; ++j) {
        X[i*n+j] = ((F[i*n+j] - (d[i]-s)*Y[i*n+j] - u[i]*Y[(i+1)*n+j]) - l[j-1]*X[i*n+(j-1)]) / w[j];
    }

    for (i = 1; i < n-1; ++i) {
        X[i*n] = (F[i*n] - l[i-1]*Y[(i-1)*n] - (d[i]-s)*Y[i*n] - u[i]*Y[(i+1)*n]) / w[0];
        for (j = 1; j < n; ++j) {
            X[i*n+j] = ((F[i*n+j] - l[i-1]*Y[(i-1)*n+j] - (d[i]-s)*Y[i*n+j] - u[i]*Y[(i+1)*n+j])
                       - l[j-1]*X[i*n+(j-1)]) / w[j];
        }
    }

    i = n-1;
    X[i*n] = (F[i*n] - l[i-1]*Y[(i-1)*n] - (d[i]-s)*Y[i*n]) / w[0];
    for (j = 1; j < n; ++j) {
        X[i*n+j] = ((F[i*n+j] - l[i-1]*Y[(i-1)*n+j] - (d[i]-s)*Y[i*n+j]) - l[j-1]*X[i*n+(j-1)]) / w[j];
    }

    /* Back substitution */

    for (j = 0; j < n; ++j) {
        for (i = n-2; i >= 0; --i) {
            X[j*n+i] = X[j*n+i] - v[i]*X[j*n+(i+1)];
        }
    }

    /* Now do the rest of the shifts (if any) */

    for (k = 1; k < nshifts; ++k) {

        s = p[k];

        /* Solve:
         *
         *     X = (T+p(j)*I) \ ( F' - ((T-p(j)*I)*X)' )
         */
        
        /* Forward pass
         *
         * Optimization: Even though Y should be transposed, pretend that
         *               it's not and remember it later. This allows us to
         *               have all arrays in the loops use the same smallest
         *               index, which makes everything more cache-friendly
         *               and thus faster.
         */

        w[0] = d[0] + s;
        Y[0] = (F[0] - (d[0]-s)*X[0] - u[0]*X[1]) / w[0];
        for (j = 1; j < n-1; ++j) {
            Y[j] = (F[j] - l[j-1]*X[j-1] - (d[j]-s)*X[j] - u[j]*X[j+1]) / w[0];
        }
        Y[n-1] = (F[n-1] - l[n-2]*X[n-2] - (d[n-1]-s)*X[n-1]) / w[0];

        for (i = 1; i < n; ++i) {
            v[i-1] = u[i-1] / w[i-1];
            w[i] = (d[i]+s) - l[i-1]*v[i-1];
            Y[i*n] = ((F[i*n] - (d[0]-s)*X[i*n] - u[0]*X[i*n+1]) - l[i-1]*Y[(i-1)*n]) / w[i];
            for (j = 1; j < n-1; ++j) {
                Y[i*n+j] = ((F[i*n+j] - l[j-1]*X[i*n+(j-1)] - (d[j]-s)*X[i*n+j] - u[j]*X[i*n+(j+1)])
                            - l[i-1]*Y[(i-1)*n+j]) / w[i];
            }
            Y[i*n+(n-1)] = ((F[i*n+(n-1)] - l[n-2]*X[i*n+(n-2)] - (d[n-1]-s)*X[i*n+(n-1)]) - l[i-1]*Y[(i-1)*n+(n-1)]) / w[i];
        }

        /* Back substitution */

        for (i = n-2; i >= 0; --i) {
            for (j = 0; j < n; ++j) {
                Y[i*n+j] = Y[i*n+j] - v[i]*Y[(i+1)*n+j];
            }
        }
        
        /* Solve:
         * 
         *    X = (T+p(j)*I) \ ( F  - ((T-p(j)*I)*X)' )
         */

        /* Forward pass
         *
         * Optimization: We don't need to compute w and v again, since we
         *               just did.
         *
         * Optimization: Remember that Y was stored transposed. This is
         *               actually good for us, since the mixed
         *               transpositions of the right-hand side here makes
         *               things right. However, in order for X to be
         *               correct, we must transpose the result. Naively,
         *               this would mean that X and the right-hand side
         *               can't share the same smallest index. But we can
         *               invert the loops entirely here---not only do i and
         *               j switch places, but the pre- and post-
         *               computations must be changed. This allows us to
         *               have all arrays in the loops use the same smallest
         *               index, which makes everything more cache-friendly
         *               and thus faster.
         */

        i = 0;
        X[i*n] = (F[i*n] - (d[i]-s)*Y[i*n] - u[i]*Y[(i+1)*n]) / w[0];
        for (j = 1; j < n; ++j) {
            X[i*n+j] = ((F[i*n+j] - (d[i]-s)*Y[i*n+j] - u[i]*Y[(i+1)*n+j]) - l[j-1]*X[i*n+(j-1)]) / w[j];
        }

        for (i = 1; i < n-1; ++i) {
            X[i*n] = (F[i*n] - l[i-1]*Y[(i-1)*n] - (d[i]-s)*Y[i*n] - u[i]*Y[(i+1)*n]) / w[0];
            for (j = 1; j < n; ++j) {
                X[i*n+j] = ((F[i*n+j] - l[i-1]*Y[(i-1)*n+j] - (d[i]-s)*Y[i*n+j] - u[i]*Y[(i+1)*n+j])
                           - l[j-1]*X[i*n+(j-1)]) / w[j];
            }
        }

        i = n-1;
        X[i*n] = (F[i*n] - l[i-1]*Y[(i-1)*n] - (d[i]-s)*Y[i*n]) / w[0];
        for (j = 1; j < n; ++j) {
            X[i*n+j] = ((F[i*n+j] - l[i-1]*Y[(i-1)*n+j] - (d[i]-s)*Y[i*n+j]) - l[j-1]*X[i*n+(j-1)]) / w[j];
        }

        /* Back substitution */

        for (j = 0; j < n; ++j) {
            for (i = n-2; i >= 0; --i) {
                X[j*n+i] = X[j*n+i] - v[i]*X[j*n+(i+1)];
            }
        }
    }

    mxFree((double *)w);
    mxFree((double *)v);
    mxFree((double *)Y);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double *d, *u, *l, *p, *F, *X;
    size_t n, nshifts;

    /* Check for the proper number of arguments */
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("FastPoisson:adi_tri:nrhs", "Five inputs required.");
    }

    if (nlhs != 1) {
        mexErrMsgIdAndTxt("FastPoisson:adi_tri:nlhs", "One output required.");
    }

    /* Check that p is a row vector */
    if (mxGetM(prhs[3]) != 1) {
        mexErrMsgIdAndTxt("FastPoisson:adi_tri:notRowVector", "Input must be a row vector.");
    }

    /* Check that d, u, l, p are column vectors */
    if (mxGetN(prhs[0]) != 1 ||
        mxGetN(prhs[1]) != 1 ||
        mxGetN(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("FastPoisson:adi_tri:notColumnVector", "Input must be a column vector.");
    }

    /* Get the dimensions */
    n = mxGetM(prhs[0]);
    nshifts = mxGetN(prhs[3]);

    /* Check that the dimensions agree */
    if (mxGetM(prhs[1]) != n-1 ||
        mxGetM(prhs[2]) != n-1 ||
        mxGetM(prhs[4]) != n) {
        mexErrMsgIdAndTxt("FastPoisson:adi_tri:mismatchedSizes", "Dimensions must agree.");
    }

    /* Check that F is a square matrix */
    if (mxGetN(prhs[4]) != n) {
        mexErrMsgIdAndTxt("FastPoisson:adi_tri:notSquareMatrix", "Input must be a square matrix.");
    }

    /* Get the input data pointers */
    d = mxGetPr(prhs[0]);
    u = mxGetPr(prhs[1]);
    l = mxGetPr(prhs[2]);
    p = mxGetPr(prhs[3]);
    F = mxGetPr(prhs[4]);

    /* Create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)n, (mwSize)n, mxREAL);

    /* Get the output data pointer */
    X = mxGetPr(plhs[0]);

    /* Do the computation */
    adi_tri(d, u, l, (mwSize)n, (mwSize)nshifts, p, F, X);
}