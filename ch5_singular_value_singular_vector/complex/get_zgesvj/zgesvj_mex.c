/* MATLAB R2025b MEX wrapper for LAPACK ZGESVJ.
 *
 * [U,S,V,sva,rwork,info] = zgesvj_mex(A,joba)
 *
 * A must be a full complex-double M-by-N matrix with M >= N. JOBU='U'
 * and JOBV='V' are fixed because the experiments always need a thin SVD.
 */

#include "mex.h"
#include "lapack.h"
#include <ctype.h>
#include <stddef.h>
#include <string.h>

static char read_joba(const mxArray *arg)
{
    char text[8] = {0};
    char joba;
    if (!mxIsChar(arg) || mxGetString(arg,text,sizeof(text)) != 0 || text[0] == '\0')
        mexErrMsgIdAndTxt("zgesvj_mex:joba","JOBA must be one of 'G', 'L', or 'U'.");
    joba = (char)toupper((unsigned char)text[0]);
    if (strchr("GLU",joba) == NULL)
        mexErrMsgIdAndTxt("zgesvj_mex:joba","JOBA must be one of 'G', 'L', or 'U'.");
    return joba;
}

static void give_output(int nlhs,mxArray *plhs[],int index,mxArray *value)
{
    if (index < nlhs) plhs[index] = value;
    else mxDestroyArray(value);
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    const mxArray *Ain;
    mxArray *Uarr,*Sarr,*Varr,*SVAarr,*RWORKarr,*INFOarr;
    mxComplexDouble *cwork;
    double *A,*V,*S,*sva,*rwork;
    mwSize mmw,nmw;
    ptrdiff_t m,n,lda,mv,ldv,lwork,lrwork,info = 0;
    char joba,jobu = 'U',jobv = 'V';
    ptrdiff_t k;

    if (nrhs != 2 || nlhs > 6)
        mexErrMsgIdAndTxt("zgesvj_mex:arity",
            "Usage: [U,S,V,sva,rwork,info] = zgesvj_mex(A,joba)");

    Ain = prhs[0];
    if (!mxIsDouble(Ain) || !mxIsComplex(Ain) || mxIsSparse(Ain))
        mexErrMsgIdAndTxt("zgesvj_mex:type","A must be a full complex-double matrix.");

    mmw = mxGetM(Ain);
    nmw = mxGetN(Ain);
    if (mmw < nmw || nmw == 0)
        mexErrMsgIdAndTxt("zgesvj_mex:shape","Require M >= N >= 1.");

    joba = read_joba(prhs[1]);
    m = (ptrdiff_t)mmw;
    n = (ptrdiff_t)nmw;
    lda = m;
    mv = 0; /* Ignored by LAPACK when JOBV='V'. */
    ldv = n;
    lwork = m+n;
    lrwork = n > 6 ? n : 6;

    Uarr = mxDuplicateArray(Ain);
    Varr = mxCreateNumericMatrix(nmw,nmw,mxDOUBLE_CLASS,mxCOMPLEX);
    SVAarr = mxCreateDoubleMatrix(nmw,1,mxREAL);
    RWORKarr = mxCreateDoubleMatrix((mwSize)lrwork,1,mxREAL);
    cwork = (mxComplexDouble *)mxCalloc((size_t)lwork,sizeof(mxComplexDouble));

    A = (double *)mxGetComplexDoubles(Uarr);
    V = (double *)mxGetComplexDoubles(Varr);
    sva = mxGetDoubles(SVAarr);
    rwork = mxGetDoubles(RWORKarr);

    zgesvj(&joba,&jobu,&jobv,&m,&n,A,&lda,sva,&mv,V,&ldv,
            (double *)cwork,&lwork,rwork,&lrwork,&info);
    mxFree(cwork);

    Sarr = mxCreateDoubleMatrix(nmw,nmw,mxREAL);
    S = mxGetDoubles(Sarr);
    for (k=0;k<n;k++) S[k+k*n] = rwork[0]*sva[k];

    INFOarr = mxCreateDoubleScalar((double)info);
    give_output(nlhs,plhs,0,Uarr);
    give_output(nlhs,plhs,1,Sarr);
    give_output(nlhs,plhs,2,Varr);
    give_output(nlhs,plhs,3,SVAarr);
    give_output(nlhs,plhs,4,RWORKarr);
    give_output(nlhs,plhs,5,INFOarr);
}
