/* MATLAB R2025b MEX wrapper for LAPACK CGEJSV.
 *
 * [U,S,V,sva,rwork,iwork,info] = cgejsv_mex(A)
 * Fixed options: JOBA='C', JOBU='U', JOBV='V', JOBR='R', JOBT='N', JOBP='N'.
 */

#include "mex.h"
#include "lapack.h"
#include <stddef.h>

static void give_output(int nlhs,mxArray *plhs[],int index,mxArray *value)
{
    if (index < nlhs) plhs[index] = value;
    else mxDestroyArray(value);
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    const mxArray *Ain;
    mxArray *Awork,*Uarr,*Sarr,*Varr,*SVAarr,*RWORKarr,*IWORKarr,*INFOarr;
    mxComplexSingle *cwork;
    ptrdiff_t *iwork;
    float *A,*U,*V,*S,*sva,*rwork;
    double *iworkOut;
    mwSize mmw,nmw;
    ptrdiff_t m,n,lda,ldu,ldv,lwork,lrwork,liwork,info = 0,k;
    char joba='C',jobu='U',jobv='V',jobr='R',jobt='N',jobp='N';

    if (nrhs != 1 || nlhs > 7)
        mexErrMsgIdAndTxt("cgejsv_mex:arity",
            "Usage: [U,S,V,sva,rwork,iwork,info] = cgejsv_mex(A)");

    Ain = prhs[0];
    if (!mxIsSingle(Ain) || !mxIsComplex(Ain) || mxIsSparse(Ain))
        mexErrMsgIdAndTxt("cgejsv_mex:type","A must be a full complex-single matrix.");
    mmw = mxGetM(Ain);
    nmw = mxGetN(Ain);
    if (mmw < nmw || nmw == 0)
        mexErrMsgIdAndTxt("cgejsv_mex:shape","Require M >= N >= 1.");

    m=(ptrdiff_t)mmw; n=(ptrdiff_t)nmw; lda=m; ldu=m; ldv=n;
    lwork = 5*n+2*n*n;
    if (lwork < m+n) lwork = m+n;
    /* Use conservative workspaces accepted by MATLAB's bundled LAPACK.
     * Some GEJSV implementations require the larger preprocessor workspace
     * even when the selected option set has a smaller documented minimum. */
    lrwork = 2*m > n ? 2*m : n;
    if (lrwork < 7) lrwork = 7;
    liwork = m+3*n;
    if (liwork < 4) liwork = 4;

    Awork = mxDuplicateArray(Ain);
    Uarr = mxCreateNumericMatrix(mmw,nmw,mxSINGLE_CLASS,mxCOMPLEX);
    Varr = mxCreateNumericMatrix(nmw,nmw,mxSINGLE_CLASS,mxCOMPLEX);
    SVAarr = mxCreateNumericMatrix(nmw,1,mxSINGLE_CLASS,mxREAL);
    RWORKarr = mxCreateNumericMatrix((mwSize)lrwork,1,mxSINGLE_CLASS,mxREAL);
    cwork = (mxComplexSingle *)mxCalloc((size_t)lwork,sizeof(mxComplexSingle));
    iwork = (ptrdiff_t *)mxCalloc((size_t)liwork,sizeof(ptrdiff_t));

    A = (float *)mxGetComplexSingles(Awork);
    U = (float *)mxGetComplexSingles(Uarr);
    V = (float *)mxGetComplexSingles(Varr);
    sva = mxGetSingles(SVAarr);
    rwork = mxGetSingles(RWORKarr);

    cgejsv(&joba,&jobu,&jobv,&jobr,&jobt,&jobp,&m,&n,A,&lda,sva,
            U,&ldu,V,&ldv,(float *)cwork,&lwork,rwork,&lrwork,iwork,&info);
    mxFree(cwork);
    mxDestroyArray(Awork);

    Sarr = mxCreateNumericMatrix(nmw,nmw,mxSINGLE_CLASS,mxREAL);
    S = mxGetSingles(Sarr);
    {
        float scale = rwork[1] != 0.0f ? rwork[0]/rwork[1] : 1.0f;
        for (k=0;k<n;k++) S[k+k*n] = scale*sva[k];
    }

    IWORKarr = mxCreateDoubleMatrix((mwSize)liwork,1,mxREAL);
    iworkOut = mxGetDoubles(IWORKarr);
    for (k=0;k<liwork;k++) iworkOut[k] = (double)iwork[k];
    mxFree(iwork);

    INFOarr = mxCreateDoubleScalar((double)info);
    give_output(nlhs,plhs,0,Uarr);
    give_output(nlhs,plhs,1,Sarr);
    give_output(nlhs,plhs,2,Varr);
    give_output(nlhs,plhs,3,SVAarr);
    give_output(nlhs,plhs,4,RWORKarr);
    give_output(nlhs,plhs,5,IWORKarr);
    give_output(nlhs,plhs,6,INFOarr);
}
