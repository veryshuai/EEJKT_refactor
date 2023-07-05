/*==========================================================
 * new_vec.c
 *
 * reads in the number of simulated firms, the number of demand shock types,
 * and the distribution of demand shock types.  The output assigns the number of
 * new clients randomly to a horizontal vector with size the number of types.  
 * 
 * The calling syntax is:
 *
 *		outVec = new_vec(number of simulated firms, number of demand shock types, cdf of demand shock types)
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include "mex.h"
#include <iostream>
using namespace std;

/* Draw a uniform rand */
float drawUniformRand(){
    
    float uniRand; 
    float randMaxSingle;
            
    /* this makes RAND MAX, an integer, into a float for dividing */
    randMaxSingle = RAND_MAX;
    
    uniRand = (((float) rand())/ randMaxSingle); 
    
    return uniRand;
}

/* Find the first instance that a number is greater than rand, and return the index */
int typeIndex(float uniRand, int Zsz, float *CdfZ){
    
    int index; 
    
    index = 1;
    for (int i=1; (uniRand > CdfZ[i-1]) & (i <= Zsz); i++){
        index = i+1;        
    }
    
    return index;
}

/* The computational routine */
void getTypeNumbers(int client_number, int Zsz, float *outMat, float *CdfZ)
{
    
    float uniRand;
    int newTypeInd;
    
    for (int c=1; c<=client_number; c++) {
        uniRand = drawUniformRand();
        newTypeInd = typeIndex(uniRand,Zsz,CdfZ);
        outMat[newTypeInd-1] = outMat[newTypeInd-1] + 1;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize client_number;              /* input scalar */
    mwSize Zsz;                        /* input scalar */
    mwSize *Zsz_array;                    /* one dimensional array */
    float *inMatrix;               /* 1xN input matrix */
    float *outMatrix;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first two arguments are scalar */
    if( !mxIsSingle(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Number of clients must be a scalar.");
    }
    
    if( !mxIsSingle(prhs[1]) || 
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Number of demand shock types must be a scalar.");
    }
    
    /* make sure the second input argument is type float */
    if( !mxIsSingle(prhs[2]) || 
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notSingle","CDF vector must be type float.");
    }

    /* get the value of the scalar inputs  */
    client_number = mxGetScalar(prhs[0]);
    Zsz = mxGetScalar(prhs[1]);
    Zsz_array[1] = {Zsz};
    
    /* check that number of rows in second input argument is 1 */
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a column vector of with size number of demand shocks.");
    }
    
    /* create a pointer to the real data in the input matrix
    #if MX_HAS_INTERLEAVED_COMPLEX
    */
    inMatrix = (float*)mxGetData(prhs[2]);
    /*
    #else
    inMatrix = mxGetPr(prhs[2]);
    #endif
    */

    /* create the output matrix */
    /* plhs[0] = mxCreateDoubleMatrix(1,Zsz,mxREAL); */
    plhs[0] = mxCreateNumericArray(1,Zsz_array,mxSINGLE_CLASS,mxREAL);

    /* get a pointer to the real data in the output matrix
        #if MX_HAS_INTERLEAVED_COMPLEX
    */
    outMatrix = (float*)mxGetData(plhs[0]);
    /* for very old versions of matlab (pre 2017)
    #else
    outMatrix = mxGetPr(plhs[0]);
    #endif
    */

    /* call the computational routine */
    getTypeNumbers(client_number,Zsz,outMatrix,inMatrix);
}
