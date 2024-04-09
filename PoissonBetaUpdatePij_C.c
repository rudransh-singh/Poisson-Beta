/* PoissonBetaUpdatePij_C.c
*
* By Jong Kyoung Kim 
* jkkim@ebi.ac.uk
* Last Update: September 24 2012
*/

#include "mex.h"
#include "PoissonBetaCparam.c"

#define COUNTDATA prhs[0]
#define CPARAMIN prhs[1]
#define APARAMIN prhs[2]
#define CPARAMOUT_PIJ plhs[0]
#define MAXSTEPSIZE 1000

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    Cparam *cparam;
    int ii, jj, J, K;
    double log_Px, xx, xprime, widths, log_uprime, rr, x_l, x_r;
    double *count_data, *size_factor;
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
    cparam->Si = mxGetPr(mxGetField(CPARAMIN,0,"Si"));
    cparam->Koni = mxGetPr(mxGetField(CPARAMIN,0,"Koni"));
    cparam->Koffi = mxGetPr(mxGetField(CPARAMIN,0,"Koffi"));
        
    /* Read Pij */
    cparam->Pij = mxMalloc(sizeof(double*)*cparam->num_gene);
    mdouble = mxGetPr(mxGetField(CPARAMIN,0,"Pij"));
    for (ii = 0; ii < cparam->num_gene; ii++) {
        cparam->Pij[ii] = mxMalloc(sizeof(double)*cparam->num_replicate);
        for (jj = 0; jj < cparam->num_replicate; jj++) {
            cparam->Pij[ii][jj] = mdouble[jj*cparam->num_gene+ii];
        }
    }          

    /* Read COUNTDATA */ 
    count_data = mxGetPr(COUNTDATA);
    
    /* Read APARAMIN */
    size_factor = mxGetPr(mxGetField(APARAMIN,0,"size_factor"));
    
    /*init_genrand64((unsigned long)time(NULL));*/
    for (ii = 0; ii < cparam->num_gene; ii++) {
        for (jj = 0; jj < cparam->num_replicate; jj++) {
            xx = cparam->Pij[ii][jj];                
            widths = xx/2.0;
            log_Px = logpbeta(cparam->Koni[ii], cparam->Koffi[ii], count_data[jj*cparam->num_gene+ii], size_factor[jj*cparam->num_gene+ii], cparam->Si[ii], xx);
            log_uprime = log(my_rand()) + log_Px;
            xprime = xx;
            rr = my_rand();
            x_l = xx - rr*widths;
            x_r = xx + (1.0-rr)*widths;
            J = (int)floor(my_rand()*MAXSTEPSIZE);
            K = MAXSTEPSIZE-1-J; 
            while ((logpbeta(cparam->Koni[ii], cparam->Koffi[ii], count_data[jj*cparam->num_gene+ii], size_factor[jj*cparam->num_gene+ii], cparam->Si[ii], x_l) > log_uprime) && (x_l - widths >= 0.0) && J > 0) {
                x_l = x_l - widths;
                J = J - 1;
            }
            while ((logpbeta(cparam->Koni[ii], cparam->Koffi[ii], count_data[jj*cparam->num_gene+ii], size_factor[jj*cparam->num_gene+ii], cparam->Si[ii], x_r) > log_uprime) && (x_r + widths <= 1.0) && K > 0) {
                x_r = x_r + widths;
                K = K - 1;
            }
            if (x_l < 0) {
                x_l = 0.0;
            }
            if (x_r > 1.0) {
                x_r = 1.0;
            }
            while (1) {
                xprime = my_rand()*(x_r - x_l) + x_l;
                log_Px = logpbeta(cparam->Koni[ii], cparam->Koffi[ii], count_data[jj*cparam->num_gene+ii], size_factor[jj*cparam->num_gene+ii], cparam->Si[ii], xprime);
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
            cparam->Pij[ii][jj] = xprime;
        }
    } 
    /* Write Pij */
    CPARAMOUT_PIJ = mxCreateDoubleMatrix(cparam->num_gene,cparam->num_replicate,mxREAL);
    mdouble = mxGetPr(CPARAMOUT_PIJ);
    for (ii = 0; ii < cparam->num_gene; ii++) {
        for (jj = 0; jj < cparam->num_replicate; jj++) {
            mdouble[jj*cparam->num_gene+ii] = cparam->Pij[ii][jj];
        }
    }    
    for (ii = 0; ii < cparam->num_gene; ii++) {
        mxFree(cparam->Pij[ii]);                
    }   
    mxFree(cparam->Pij);
    mxFree(cparam);
}
