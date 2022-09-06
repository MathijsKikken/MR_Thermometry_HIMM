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
 * Outputs:
 * outMatrix - the result of the operation above
 * change - the squared change between outMatrix and inMatrix
 *
 * The calling syntax is:
 *
 * [outMatrix, change] = MeanAdjacentWithMaskIncChange(inMatrix, mask)
 *
 *========================================================*/

#include "mex.h"

/* The computational routine */
void MeanAdjacentWithMask(double *inMatrix, bool *inMask, double *outMatrix, double *change, mwSize n, mwSize m)
{
    mwSize i;
    int count;
    double tmp;
    
    // Loop over all
    for (int i=0; i<n*m; i++) {
        if (inMask[i] == 1) {
            // If on boundary or in mask, don't change
            outMatrix[i] = inMatrix[i];
        } else {
            // Else do
            // First check which nearby you can use
            count = 0;
            tmp = 0;
            if (i >= m) {
                count++;
                tmp += inMatrix[i-m];
            }
            if (i < n*m - m) {
                count++;
                tmp += inMatrix[i+m];
            }
            if (i % m > 0) {
                count++;
                tmp += inMatrix[i-1];
            }
            if (i % m < m - 1) {
                count++;
                tmp += inMatrix[i+1];
            }
            // Calculate new value
            outMatrix[i] = tmp / count;
            // Also calculate squared change
            *change += (inMatrix[i] - outMatrix[i]) * (inMatrix[i] - outMatrix[i]);
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
    size_t nrows;                   /* size of matrix */
    double *outMatrix;              /* output matrix */
    double *change;                 /* change between input and output (squared)*/

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two outputs required.");
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
    
    /* check that number of rows in second input argument equals number of columns */
    if(mxGetM(prhs[1])!= mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input 2 must be a matrix.");
    }
    
    /* check that both matrices have same size */
    if(mxGetN(prhs[1])!= mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Both matrices must be same size.");
    }

    /* create a pointer to the real data in the input matrix  */
    #if MX_HAS_INTERLEAVED_COMPLEX
    inMatrix = mxGetDoubles(prhs[0]);
    #else
    inMatrix = mxGetPr(prhs[0]);
    #endif
    
    /* create a pointer to the real data in the mask matrix  */
    inMask = mxGetLogicals(prhs[1]);
    

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[0]);
    nrows = mxGetM(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nrows,(mwSize)ncols,mxREAL);
    /* create the secondary output: change */
    plhs[1] = mxCreateDoubleScalar(0.0);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    outMatrix = mxGetDoubles(plhs[0]);
    #else
    outMatrix = mxGetPr(plhs[0]);
    #endif
  
    /* get a pointer to the real data of the second result: change */
    #if MX_HAS_INTERLEAVED_COMPLEX
    change = mxGetDoubles(plhs[1]);
    #else
    change = mxGetPr(plhs[1]);
    #endif
        
    *change = 0; // Predefine change as 0
        
    /* call the computational routine */
    MeanAdjacentWithMask(inMatrix, inMask, outMatrix, change, (mwSize)ncols, (mwSize)nrows);
}
