/*==========================================================
 * MeanAdjacentWithMaskIncChange.c
 *
 * Inputs:
 * inMatrix - the input matrix
 * inMask - a mask that determines if that location in the matrix can be changes
 *
 * Task:
 * Changes the value of each location in the matrix to the mean of its four neighbors
 * Unless it is on the boundary or in the mask
 *
 * Output:
 * outMatrix - the result of the operation above
 *
 * The calling syntax is:
 *
 * outMatrix = MeanAdjacentWithMaskIncChange(inMatrix, mask)
 *
 *========================================================*/

#include "mex.h"

/* The computational routine */
void MeanAdjacentWithMask(double *inMatrix, bool *inMask, double *outMatrix, mwSize n)
{
    mwSize i;
    // Loop over all
    for (int i=0; i<n*n; i++) {
        if (i < n || i > n*n - n || i % n == 0 || i % n == n - 1 || inMask[i] == 1) {
            // If on boundary or in mask, don't change
            outMatrix[i] = inMatrix[i];
        } else {
            // Else do
            outMatrix[i] = (inMatrix[i-1] + inMatrix[i+1] + inMatrix[i-n] + inMatrix[i+n])/4;
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *inMatrix;               /* NxN input matrix */
    bool *inMask;                   /* NxN mask matrix */
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    
    /* make sure the first input argument is type double */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix 1 must be double.");
    }
    
    /* make sure the second input argument is type logical */
    if( !mxIsLogical(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix 2 must be logical.");
    }
    
    /* check that number of rows in first input argument equals number of columns */
    if(mxGetM(prhs[0])!= mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input 2 must be a matrix.");
    }
    
    /* check that number of rows in second input argument equals number of columns */
    if(mxGetM(prhs[1])!= mxGetN(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input 2 must be a matrix.");
    }
    
    /* check that both matrices have same size */
    if(mxGetM(prhs[1])!= mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Both matrices must be same size.");
    }

    /* create a pointer to the real data in the input matrix  */
    #if MX_HAS_INTERLEAVED_COMPLEX
    inMatrix = mxGetDoubles(prhs[0]);
    #else
    inMatrix = mxGetPr(prhs[0]);
    #endif
    
    /* create a pointer to the logical data in the mask matrix  */
    inMask = mxGetLogicals(prhs[1]);

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)ncols,(mwSize)ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    outMatrix = mxGetDoubles(plhs[0]);
    #else
    outMatrix = mxGetPr(plhs[0]);
    #endif

    /* call the computational routine */
    MeanAdjacentWithMask(inMatrix, inMask, outMatrix,(mwSize)ncols);
}
