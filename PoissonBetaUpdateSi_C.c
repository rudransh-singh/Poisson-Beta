/* PoissonBetaUpdateSi_C.c
*
* By Jong Kyoung Kim 
* jkkim@ebi.ac.uk
* Last Update: July 10 2012
*/

#include "mex.h"
#include "PoissonBetaCparam.c"

#define CPARAMIN prhs[2]
#define APARAMIN prhs[1]
#define SUM_COUNT prhs[0]
#define CPARAMOUT_SI plhs[0]
#define MAXSTEPSIZE 1000

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    Cparam *cparam;
    int ii, jj, J, K;
    double log_Px, xx, xprime, widths, log_uprime, rr, x_l, x_r;
    double *alpha_s, *beta_s, *size_factor, *sum_count, *Pij;
    double *mdouble;
        
    /* Check for proper number of arguments. */
    if (nrhs != 3) {
      mexErrMsgTxt("3 inputs required.");
    } else if (nlhs != 1) {
      mexErrMsgTxt("1 output required");      
    }
    
    /* Read CPARAMIN */ 
    cparam = mxMalloc(sizeof(Cparam));
    
    cparam->num_gene = (int)mxGetScalar(mxGetField(CPARAMIN,0,"num_gene"));
    cparam->num_replicate = (int)mxGetScalar(mxGetField(CPARAMIN,0,"num_replicate"));     
        
    /* Read Si */
    cparam->Si = mxMalloc(sizeof(double)*cparam->num_gene);
    mdouble = mxGetPr(mxGetField(CPARAMIN,0,"Si"));
    for (ii = 0; ii < cparam->num_gene; ii++) {
        cparam->Si[ii] = mdouble[ii];  
    }         
    
    /* Read Pij */
    Pij = mxGetPr(mxGetField(CPARAMIN,0,"Pij"));
    
    /* Read sum_count */
    sum_count = mxGetPr(SUM_COUNT);
    
    /* Read APARAMIN */
    size_factor = mxGetPr(mxGetField(APARAMIN,0,"size_factor"));
    alpha_s = mxGetPr(mxGetField(APARAMIN,0,"alpha_si"));
    beta_s = mxGetPr(mxGetField(APARAMIN,0,"beta_si"));
    
    for (ii = 0; ii < cparam->num_gene; ii++) {
        xx = cparam->Si[ii];
        widths = cparam->Si[ii]/2.0;
        log_Px = logpgammas(alpha_s[ii], beta_s[ii], size_factor, Pij, sum_count[ii], cparam->num_gene, cparam->num_replicate, ii, xx);
        log_uprime = log(my_rand()) + log_Px;
        xprime = xx;
        rr = my_rand();
        x_l = xx - rr*widths;
        x_r = xx + (1-rr)*widths;
        J = (int)floor(my_rand()*MAXSTEPSIZE);
        K = MAXSTEPSIZE-1-J;
        while (logpgammas(alpha_s[ii], beta_s[ii], size_factor, Pij, sum_count[ii], cparam->num_gene, cparam->num_replicate, ii, x_l) > log_uprime && x_l - widths >= 0.0 && J > 0) {
            x_l = x_l - widths;
            J = J - 1;
        }
        while (logpgammas(alpha_s[ii], beta_s[ii], size_factor, Pij, sum_count[ii], cparam->num_gene, cparam->num_replicate, ii, x_r) > log_uprime && K > 0) {
            x_r = x_r + widths;
            K = K - 1;
        }
        while (1) {
            xprime = my_rand()*(x_r - x_l) + x_l;
            log_Px = logpgammas(alpha_s[ii], beta_s[ii], size_factor, Pij, sum_count[ii], cparam->num_gene, cparam->num_replicate, ii, xprime);
            if (log_Px > log_uprime) {
                break;
            }
            else {
                if (xprime > xx) {
                    x_r = xprime;
                }
                else if (xprime < xx) {
                    x_l = xprime;
                }
                else {
                    mexErrMsgTxt("BUG DETECTED: Shrunk to current position and still not acceptable.");
                }
            }
        }
        cparam->Si[ii] = xprime;
    } 
    /* Write Si */
    CPARAMOUT_SI = mxCreateDoubleMatrix(cparam->num_gene, 1, mxREAL);
    mdouble = mxGetPr(CPARAMOUT_SI);
    for (ii = 0; ii < cparam->num_gene; ii++) {
        mdouble[ii] = cparam->Si[ii];
    }
    mxFree(cparam->Si);
    
    mxFree(cparam);
}