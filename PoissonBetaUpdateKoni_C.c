/* PoissonBetaUpdateKoni_C.c
*
* By Jong Kyoung Kim 
* jkkim@ebi.ac.uk
* Last Update: July 10 2012
*/

#include "mex.h"
#include "PoissonBetaCparam.c"

#define CPARAMIN prhs[1]
#define APARAMIN prhs[0]
#define CPARAMOUT_KONI plhs[0]
#define MAXSTEPSIZE 1000

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    Cparam *cparam;
    int ii, jj, J, K;
    double log_Px, xx, xprime, widths, log_uprime, rr, x_l, x_r;
    double *alpha_kon, *beta_kon;
    double *mdouble, *lnoneminusPij;
        
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
    
    cparam->Koffi = mxGetPr(mxGetField(CPARAMIN,0,"Koffi"));
        
    /* Read Koni */
    cparam->Koni = mxMalloc(sizeof(double)*cparam->num_gene);
    mdouble = mxGetPr(mxGetField(CPARAMIN,0,"Koni"));
    for (ii = 0; ii < cparam->num_gene; ii++) {
        cparam->Koni[ii] = mdouble[ii];  
    }         
    
    /* Read lnoneminusPij */
    lnoneminusPij = mxGetPr(mxGetField(CPARAMIN,0,"lnoneminusPij"));

    /* Read APARAMIN */
    alpha_kon = mxGetPr(mxGetField(APARAMIN,0,"alpha_koni"));
    beta_kon = mxGetPr(mxGetField(APARAMIN,0,"beta_koni"));
    
    for (ii = 0; ii < cparam->num_gene; ii++) {
        xx = cparam->Koni[ii];
        widths = cparam->Koni[ii]/2.0;
        log_Px = logpgammakon(alpha_kon[ii], beta_kon[ii], cparam->Koffi[ii], lnoneminusPij, cparam->num_gene, cparam->num_replicate, ii, xx);
        log_uprime = log(my_rand()) + log_Px;
        xprime = xx;
        rr = my_rand();
        x_l = xx - rr*widths;
        x_r = xx + (1-rr)*widths;
        J = (int)floor(my_rand()*MAXSTEPSIZE);
        K = MAXSTEPSIZE-1-J;
        while (logpgammakon(alpha_kon[ii], beta_kon[ii], cparam->Koffi[ii], lnoneminusPij, cparam->num_gene, cparam->num_replicate, ii, x_l) > log_uprime && x_l - widths >= 0.0 && J > 0) {
            x_l = x_l - widths;
            J = J - 1;
        }
        while (logpgammakon(alpha_kon[ii], beta_kon[ii], cparam->Koffi[ii], lnoneminusPij, cparam->num_gene, cparam->num_replicate, ii, x_r) > log_uprime && K > 0) {
            x_r = x_r + widths;
            K = K - 1;
        }
        while (1) {
            xprime = my_rand()*(x_r - x_l) + x_l;
            log_Px = logpgammakon(alpha_kon[ii], beta_kon[ii], cparam->Koffi[ii], lnoneminusPij, cparam->num_gene, cparam->num_replicate, ii, xprime);
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
        cparam->Koni[ii] = xprime;
    } 
    /* Write Koni */
    CPARAMOUT_KONI = mxCreateDoubleMatrix(cparam->num_gene, 1, mxREAL);
    mdouble = mxGetPr(CPARAMOUT_KONI);
    for (ii = 0; ii < cparam->num_gene; ii++) {
        mdouble[ii] = cparam->Koni[ii];
    }
    mxFree(cparam->Koni);
    
    mxFree(cparam);
}