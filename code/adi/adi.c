/*==========================================================
 * adi.c
 *
 * Solves the equation AX - XB' = F, where A and B are special
 * pentadiagonal matrices with three nonzero diagonals.
 *
 * Inputs:
 *     da (db) = main diagonal of A (B)
 *     ua (ub) = upper diagonal of A (B)
 *     la (lb) = lower diagonal of A (B)
 *     p, q = shift parameters
 *     F = right-hand side matrix
 *
 * Outputs:
 *     X, where X is the result of
 *
 *     X = zeros(n,n);
 *     for j = 1:numel(p)
 *         X = (B+q(j)*In) \ ( F' - ((A+q(j)*Im)*X)' );
 *         X = (A+p(j)*Im) \ ( F - ((B+p(j)*In)*X)' );
 *     end
 *
 * The calling syntax is:
 *
 *     X = adi( da, ua, la, db, ub, lb, p, q, F )
 *
 *========================================================*/

#include "mex.h"

void adi(double *da, double *ua, double *la,
         double *db, double *ub, double *lb,
         mwSize m, mwSize n, double *p, double *q, mwSize nshifts, double *F, double *X)
{
    int i, j, k;
    double s, *wa, *va, *wb, *vb, *Y;

    wa = (double *)mxMalloc(m * sizeof(double));
    va = (double *)mxMalloc(m * sizeof(double));
    wb = (double *)mxMalloc(n * sizeof(double));
    vb = (double *)mxMalloc(n * sizeof(double));
    Y  = (double *)mxMalloc(m * n * sizeof(double));

    /* Optimization: Assuming X = 0 initially lets us simplify the first
     *               iteration without explicitly setting X to zero.
     *
     * See main loop for notes on the optimizations performed here.
     */

    /* Solve:
     *  
     *     X = (B+q(j)*In) \ F'
     */

    s = q[0];

    Y[0]   = F[0]   / (db[0]+s);
    Y[m]   = F[m]   / (db[1]+s);
    Y[1]   = F[1]   / (db[0]+s);
    Y[m+1] = F[m+1] / (db[1]+s);
    for (j = 2; j < m-2; ++j) {
        Y[j]   = F[j]   / (db[0]+s);
        Y[m+j] = F[m+j] / (db[1]+s);
    }
    Y[m-2]   = F[m-2]   / (db[0]+s);
    Y[2*m-2] = F[2*m-2] / (db[1]+s);
    Y[m-1]   = F[m-1]   / (db[0]+s);
    Y[2*m-1] = F[2*m-1] / (db[1]+s);

    vb[0] = ub[0] / (db[0]+s);
    vb[1] = ub[1] / (db[1]+s);
    for (i = 2; i < n-2; ++i) {
        wb[i-2] = (db[i]+s) - lb[i-2] * vb[i-2];
        vb[i] = ub[i] / wb[i-2];
        Y[i*m]   = (F[i*m]   - lb[i-2]*Y[(i-2)*m])   / wb[i-2];
        Y[i*m+1] = (F[i*m+1] - lb[i-2]*Y[(i-2)*m+1]) / wb[i-2];
        for (j = 2; j < m-2; ++j) {
            Y[i*m+j] = (F[i*m+j] - lb[i-2]*Y[(i-2)*m+j]) / wb[i-2];
        }
        Y[i*m+(m-2)] = (F[i*m+(m-2)] - lb[i-2]*Y[(i-2)*m+(m-2)]) / wb[i-2];
        Y[i*m+(m-1)] = (F[i*m+(m-1)] - lb[i-2]*Y[(i-2)*m+(m-1)]) / wb[i-2];
    }

    for (i = n-2; i < n; ++i) {
        wb[i-2] = (db[i]+s) - lb[i-2] * vb[i-2];
        Y[i*m]   = (F[i*m]   - lb[i-2]*Y[(i-2)*m])   / wb[i-2];
        Y[i*m+1] = (F[i*m+1] - lb[i-2]*Y[(i-2)*m+1]) / wb[i-2];
        for (j = 2; j < m-2; ++j) {
            Y[i*m+j] = (F[i*m+j] - lb[i-2]*Y[(i-2)*m+j]) / wb[i-2];
        }
        Y[i*m+(m-2)] = (F[i*m+(m-2)] - lb[i-2]*Y[(i-2)*m+(m-2)]) / wb[i-2];
        Y[i*m+(m-1)] = (F[i*m+(m-1)] - lb[i-2]*Y[(i-2)*m+(m-1)]) / wb[i-2];
    }

    /* Back substitution */

    for (i = n-3; i >= 0; --i) {
        for (j = 0; j < m; ++j) {
            Y[i*m+j] = Y[i*m+j] - vb[i] * Y[(i+2)*m+j];
        }
    }

    /* Solve:
     * 
     *    X = (A+p(j)*Im) \ ( F - ((B+p(j)*In)*X)' );
     */

    /* Forward pass */

    s = p[0];

    va[0] = ua[0] / (da[0]+s);
    va[1] = ua[1] / (da[1]+s);
    for (j = 2; j < m-2; ++j) {
        wa[j-2] = (da[j]+s) - la[j-2] * va[j-2];
        va[j] = ua[j] / wa[j-2];
    }
    for (j = m-2; j < m; ++j) {
        wa[j-2] = (da[j]+s) - la[j-2] * va[j-2];
    }

    for (i = 0; i < 2; ++i) {
        X[i*m]   = (F[i*m]   - (db[i]+s)*Y[i*m]   - ub[i]*Y[(i+2)*m])   / (da[0]+s);
        X[i*m+1] = (F[i*m+1] - (db[i]+s)*Y[i*m+1] - ub[i]*Y[(i+2)*m+1]) / (da[1]+s);
        for (j = 2; j < m; ++j) {
            X[i*m+j] = ((F[i*m+j] - (db[i]+s)*Y[i*m+j] - ub[i]*Y[(i+2)*m+j]) - la[j-2]*X[i*m+j-2]) / wa[j-2];
        }
    }

    for (i = 2; i < n-2; ++i) {
        X[i*m]   = (F[i*m]   - lb[i-2]*Y[(i-2)*m]   - (db[i]+s)*Y[i*m]   - ub[i]*Y[(i+2)*m])   / (da[0]+s);
        X[i*m+1] = (F[i*m+1] - lb[i-2]*Y[(i-2)*m+1] - (db[i]+s)*Y[i*m+1] - ub[i]*Y[(i+2)*m+1]) / (da[1]+s);
        for (j = 2; j < m; ++j) {
            X[i*m+j] = ((F[i*m+j] - lb[i-2]*Y[(i-2)*m+j] - (db[i]+s)*Y[i*m+j] - ub[i]*Y[(i+2)*m+j])
                       - la[j-2]*X[i*m+j-2]) / wa[j-2];
        }
    }

    for (i = n-2; i < n; ++i) {
        X[i*m]   = (F[i*m]   - lb[i-2]*Y[(i-2)*m]   - (db[i]+s)*Y[i*m])   / (da[0]+s);
        X[i*m+1] = (F[i*m+1] - lb[i-2]*Y[(i-2)*m+1] - (db[i]+s)*Y[i*m+1]) / (da[1]+s);
        for (j = 2; j < m; ++j) {
            X[i*m+j] = ((F[i*m+j] - lb[i-2]*Y[(i-2)*m+j] - (db[i]+s)*Y[i*m+j]) - la[j-2]*X[i*m+j-2]) / wa[j-2];
        }
    }

    /* Back substitution */

    for (i = 0; i < n; ++i) {
        for (j = m-3; j >= 0; --j) {
            X[i*m+j] = X[i*m+j] - va[j] * X[i*m+j+2];
        }
    }

    /* Now do the rest of the shifts (if any) */

    for (k = 0; k < nshifts; ++k) {

        s = q[k];

        /* Solve:
         *
         *     X = (B+q(j)*In) \ ( F' - ((A+q(j)*Im)*X)' );
         */

        /* Forward pass
         *
         * Optimization: Even though Y should be transposed, pretend that
         *               it's not and remember it later. This allows us to
         *               have all arrays in the loops use the same smallest
         *               index, which makes everything more cache-friendly
         *               and thus faster.
         */

        Y[0]   = (F[0]   - (da[0]+s)*X[0]   - ua[0]*X[2])   / (db[0]+s);
        Y[m]   = (F[m]   - (da[0]+s)*X[m]   - ua[0]*X[m+2]) / (db[1]+s);
        Y[1]   = (F[1]   - (da[1]+s)*X[1]   - ua[1]*X[3])   / (db[0]+s);
        Y[m+1] = (F[m+1] - (da[1]+s)*X[m+1] - ua[1]*X[m+3]) / (db[1]+s);
        for (j = 2; j < m-2; ++j) {
            Y[j]   = (F[j]   - la[j-2]*X[j-2]   - (da[j]+s)*X[j]   - ua[j]*X[j+2])   / (db[0]+s);
            Y[m+j] = (F[m+j] - la[j-2]*X[m+j-2] - (da[j]+s)*X[m+j] - ua[j]*X[m+j+2]) / (db[1]+s);
        }
        Y[m-2]   = (F[m-2]   - la[m-4]*X[m-4]   - (da[m-2]+s)*X[m-2])   / (db[0]+s);
        Y[2*m-2] = (F[2*m-2] - la[m-4]*X[2*m-4] - (da[m-2]+s)*X[2*m-2]) / (db[1]+s);
        Y[m-1]   = (F[m-1]   - la[m-3]*X[m-3]   - (da[m-1]+s)*X[m-1])   / (db[0]+s);
        Y[2*m-1] = (F[2*m-1] - la[m-3]*X[2*m-3] - (da[m-1]+s)*X[2*m-1]) / (db[1]+s);

        vb[0] = ub[0] / (db[0]+s);
        vb[1] = ub[1] / (db[1]+s);
        for (i = 2; i < n-2; ++i) {
            wb[i-2] = (db[i]+s) - lb[i-2] * vb[i-2];
            vb[i] = ub[i] / wb[i-2];
            Y[i*m]   = ((F[i*m]   - (da[0]+s)*X[i*m]   - ua[0]*X[i*m+2]) - lb[i-2]*Y[(i-2)*m])   / wb[i-2];
            Y[i*m+1] = ((F[i*m+1] - (da[1]+s)*X[i*m+1] - ua[1]*X[i*m+3]) - lb[i-2]*Y[(i-2)*m+1]) / wb[i-2];
            for (j = 2; j < m-2; ++j) {
                Y[i*m+j] = ((F[i*m+j] - la[j-2]*X[i*m+j-2] - (da[j]+s)*X[i*m+j] - ua[j]*X[i*m+j+2])
                            - lb[i-2]*Y[(i-2)*m+j]) / wb[i-2];
            }
            Y[i*m+(m-2)] = ((F[i*m+(m-2)] - la[m-4]*X[i*m+(m-4)] - (da[m-2]+s)*X[i*m+(m-2)]) - lb[i-2]*Y[(i-2)*m+(m-2)]) / wb[i-2];
            Y[i*m+(m-1)] = ((F[i*m+(m-1)] - la[m-3]*X[i*m+(m-3)] - (da[m-1]+s)*X[i*m+(m-1)]) - lb[i-2]*Y[(i-2)*m+(m-1)]) / wb[i-2];
        }

        for (i = n-2; i < n; ++i) {
            wb[i-2] = (db[i]+s) - lb[i-2] * vb[i-2];
            Y[i*m]   = ((F[i*m]   - (da[0]+s)*X[i*m]   - ua[0]*X[i*m+2]) - lb[i-2]*Y[(i-2)*m])   / wb[i-2];
            Y[i*m+1] = ((F[i*m+1] - (da[1]+s)*X[i*m+1] - ua[1]*X[i*m+3]) - lb[i-2]*Y[(i-2)*m+1]) / wb[i-2];
            for (j = 2; j < m-2; ++j) {
                Y[i*m+j] = ((F[i*m+j] - la[j-2]*X[i*m+j-2] - (da[j]+s)*X[i*m+j] - ua[j]*X[i*m+j+2])
                            - lb[i-2]*Y[(i-2)*m+j]) / wb[i-2];
            }
            Y[i*m+(m-2)] = ((F[i*m+(m-2)] - la[m-4]*X[i*m+(m-4)] - (da[m-2]+s)*X[i*m+(m-2)]) - lb[i-2]*Y[(i-2)*m+(m-2)]) / wb[i-2];
            Y[i*m+(m-1)] = ((F[i*m+(m-1)] - la[m-3]*X[i*m+(m-3)] - (da[m-1]+s)*X[i*m+(m-1)]) - lb[i-2]*Y[(i-2)*m+(m-1)]) / wb[i-2];
        }
 
        /* Back substitution */

        for (i = n-3; i >= 0; --i) {
            for (j = 0; j < m; ++j) {
                Y[i*m+j] = Y[i*m+j] - vb[i] * Y[(i+2)*m+j];
            }
        }

        /* Solve:
         * 
         *    X = (A+p(j)*Im) \ ( F - ((B+p(j)*In)*X)' );
         */

        /* Forward pass
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

        s = p[k];

        va[0] = ua[0] / (da[0]+s);
        va[1] = ua[1] / (da[1]+s);
        for (j = 2; j < m-2; ++j) {
            wa[j-2] = (da[j]+s) - la[j-2] * va[j-2];
            va[j] = ua[j] / wa[j-2];
        }
        for (j = m-2; j < m; ++j) {
            wa[j-2] = (da[j]+s) - la[j-2] * va[j-2];
        }
        
        for (i = 0; i < 2; ++i) {
            X[i*m]   = (F[i*m]   - (db[i]+s)*Y[i*m]   - ub[i]*Y[(i+2)*m])   / (da[0]+s);
            X[i*m+1] = (F[i*m+1] - (db[i]+s)*Y[i*m+1] - ub[i]*Y[(i+2)*m+1]) / (da[1]+s);
            for (j = 2; j < m; ++j) {
                X[i*m+j] = ((F[i*m+j] - (db[i]+s)*Y[i*m+j] - ub[i]*Y[(i+2)*m+j]) - la[j-2]*X[i*m+j-2]) / wa[j-2];
            }
        }

        for (i = 2; i < n-2; ++i) {
            X[i*m]   = (F[i*m]   - lb[i-2]*Y[(i-2)*m]   - (db[i]+s)*Y[i*m]   - ub[i]*Y[(i+2)*m])   / (da[0]+s);
            X[i*m+1] = (F[i*m+1] - lb[i-2]*Y[(i-2)*m+1] - (db[i]+s)*Y[i*m+1] - ub[i]*Y[(i+2)*m+1]) / (da[1]+s);
            for (j = 2; j < m; ++j) {
                X[i*m+j] = ((F[i*m+j] - lb[i-2]*Y[(i-2)*m+j] - (db[i]+s)*Y[i*m+j] - ub[i]*Y[(i+2)*m+j])
                           - la[j-2]*X[i*m+j-2]) / wa[j-2];
            }
        }

        for (i = n-2; i < n; ++i) {
            X[i*m]   = (F[i*m]   - lb[i-2]*Y[(i-2)*m]   - (db[i]+s)*Y[i*m])   / (da[0]+s);
            X[i*m+1] = (F[i*m+1] - lb[i-2]*Y[(i-2)*m+1] - (db[i]+s)*Y[i*m+1]) / (da[1]+s);
            for (j = 2; j < m; ++j) {
                X[i*m+j] = ((F[i*m+j] - lb[i-2]*Y[(i-2)*m+j] - (db[i]+s)*Y[i*m+j]) - la[j-2]*X[i*m+j-2]) / wa[j-2];
            }
        }

        /* Back substitution */

        for (i = 0; i < n; ++i) {
            for (j = m-3; j >= 0; --j) {
                X[i*m+j] = X[i*m+j] - va[j] * X[i*m+j+2];
            }
        }
    }

    mxFree((double *)wa);
    mxFree((double *)va);
    mxFree((double *)wb);
    mxFree((double *)vb);
    mxFree((double *)Y);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double *da, *ua, *la; // diagonals of A
    double *db, *ub, *lb; // diagonals of B
    double *p, *q;        // ADI shifts
    double *F, *X;
    mwSize m, n, np, nq;

    /* Check for the proper number of arguments */
    if (nrhs != 9) {
        mexErrMsgIdAndTxt("POISSON:ADI", "Five inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("POISSON:ADI", "One output required.");
    }

    /* Check that p, q are row vectors */
    if (mxGetM(prhs[6]) != 1 ||
        mxGetM(prhs[7]) != 1) {
        mexErrMsgIdAndTxt("POISSON:ADI", "Shifts must be row vectors.");
    }

    /* Check that da, ua, la, & db, ub, lb are column vectors */
    if (mxGetN(prhs[0]) != 1 ||
        mxGetN(prhs[1]) != 1 ||
        mxGetN(prhs[2]) != 1 ||
        mxGetN(prhs[3]) != 1 ||
        mxGetN(prhs[4]) != 1 ||
        mxGetN(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("POISSON:ADI", "Diagonals must be column vectors.");
    }

    /* Get the dimensions */
    m = mxGetM(prhs[0]); // da
    n = mxGetM(prhs[3]); // db
    np = mxGetN(prhs[6]);
    nq = mxGetN(prhs[7]);

    if (np != nq) {
        mexErrMsgIdAndTxt("POISSON:ADI", "Number of shifts must be equal.");
    }

    /* Check that the dimensions agree */
    if (mxGetM(prhs[1]) != m-2 || // ua
        mxGetM(prhs[2]) != m-2 || // la
        mxGetM(prhs[4]) != n-2 || // ub
        mxGetM(prhs[5]) != n-2 || // lb
        mxGetM(prhs[8]) != m   || // F
        mxGetN(prhs[8]) != n) {   // F
        mexErrMsgIdAndTxt("POISSON:ADI", "Dimensions must agree.");
    }

    /* Get the input data pointers */
    da = mxGetPr(prhs[0]);
    ua = mxGetPr(prhs[1]);
    la = mxGetPr(prhs[2]);
    db = mxGetPr(prhs[3]);
    ub = mxGetPr(prhs[4]);
    lb = mxGetPr(prhs[5]);
    p = mxGetPr(prhs[6]);
    q = mxGetPr(prhs[7]);
    F = mxGetPr(prhs[8]);

    /* Create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);

    /* Get the output data pointer */
    X = mxGetPr(plhs[0]);

    /* Do the computation */
    adi(da, ua, la, db, ub, lb, m, n, p, q, np, F, X);
}
