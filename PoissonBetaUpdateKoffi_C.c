/* PoissonBetaUpdateKoffi_C.c
*
* By Jong Kyoung Kim 
* jkkim@ebi.ac.uk
* Last Update: July 10 2012
*/

#include "mex.h"
#include "PoissonBetaCparam.c"

#define CPARAMIN prhs[1]
#define APARAMIN prhs[0]
#define CPARAMOUT_KOFFI plhs[0]
#define MAXSTEPSIZE 1000

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    Cparam *cparam;
    int ii, jj, J, K;
    double log_Px, xx, xprime, widths, log_uprime, rr, x_l, x_r;
    double *alpha_koff, *beta_koff, *lnPij;
    double *mdouble;
        
    /* Check for proper number of arguments. */
    if (nrhs != 2) {
      mexErrMsgTxt("2 inputs required.");
    } else if (nlhs != 1) {
      mexErrMsgTxt("1 output required");      
    }
    
    /* Read CPARAMIN */ 
    cparam = mxMalloc(sizeof(Cparam));
    
    cparam->num_gene = (int)mxGetScalar(mxGetField(CPARAMIN,0,"num_gene"));
    cparam->num_replicate = (int)mxGetScalar(mxGetField(CPARAMIN,0,"num_replicate"));     
    
    cparam->Koni = mxGetPr(mxGetField(CPARAMIN,0,"Koni"));
        
    /* Read Koffi */
    cparam->Koffi = mxMalloc(sizeof(double)*cparam->num_gene);
    mdouble = mxGetPr(mxGetField(CPARAMIN,0,"Koffi"));
    for (ii = 0; ii < cparam->num_gene; ii++) {
        cparam->Koffi[ii] = mdouble[ii];  
    }         
    
    /* Read lnPij */
    lnPij = mxGetPr(mxGetField(CPARAMIN,0,"lnPij"));
    
    /* Read APARAMIN */
    alpha_koff = mxGetPr(mxGetField(APARAMIN,0,"alpha_koffi"));
    beta_koff = mxGetPr(mxGetField(APARAMIN,0,"beta_koffi"));
    
    for (ii = 0; ii < cparam->num_gene; ii++) {
        xx = cparam->Koffi[ii];
        widths = cparam->Koffi[ii]/2.0;
        log_Px = logpgammakoff(alpha_koff[ii], beta_koff[ii], cparam->Koni[ii], lnPij, cparam->num_gene, cparam->num_replicate, ii, xx);
        log_uprime = log(my_rand()) + log_Px;
        xprime = xx;
        rr = my_rand();
        x_l = xx - rr*widths;
        x_r = xx + (1-rr)*widths;
        J = (int)floor(my_rand()*MAXSTEPSIZE);
        K = MAXSTEPSIZE-1-J;
        while (logpgammakoff(alpha_koff[ii], beta_koff[ii], cparam->Koni[ii], lnPij, cparam->num_gene, cparam->num_replicate, ii, x_l) > log_uprime && x_l - widths >= 0.0 && J > 0) {
            x_l = x_l - widths;
            J = J - 1;
        }
        while (logpgammakoff(alpha_koff[ii], beta_koff[ii], cparam->Koni[ii], lnPij, cparam->num_gene, cparam->num_replicate, ii, x_r) > log_uprime && K > 0) {
            x_r = x_r + widths;
            K = K - 1;
        }
        while (1) {
            xprime = my_rand()*(x_r - x_l) + x_l;
            log_Px = logpgammakoff(alpha_koff[ii], beta_koff[ii], cparam->Koni[ii], lnPij, cparam->num_gene, cparam->num_replicate, ii, xprime);
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
                    break;
                    /*mexErrMsgTxt("BUG DETECTED: Shrunk to current position and still not acceptable.");*/
                }
            }
        }
        cparam->Koffi[ii] = xprime;
    } 
    /* Write Koffi */
    CPARAMOUT_KOFFI = mxCreateDoubleMatrix(cparam->num_gene, 1, mxREAL);
    mdouble = mxGetPr(CPARAMOUT_KOFFI);
    for (ii = 0; ii < cparam->num_gene; ii++) {
        mdouble[ii] = cparam->Koffi[ii];
    }
    mxFree(cparam->Koffi);
    
    mxFree(cparam);
}