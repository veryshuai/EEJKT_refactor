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
double drawUniformRand(){
    
    double uniRand; 
    double randMaxDouble;
            
    /* this makes RAND MAX, an integer, into a double for dividing */
    randMaxDouble = RAND_MAX;
    
    uniRand = (((double) rand())/ randMaxDouble); 
    
    return uniRand;
} 

/* Find the first instance that a number is greater than rand, and return the index */
int typeIndex(double uniRand, int Zsz, double *CdfZ){
    
    int index; 
    
    index = 1;
    for (int i=1; (uniRand > CdfZ[i-1]) & (i <= Zsz); i++){
        index = i+1;        
    }
    
    return index;
}

/* The computational routine */
void getTypeNumbers(int client_number, int Zsz, double *outMat, double *CdfZ)
{
    
    double uniRand;
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
    int client_number;              /* input scalar */
    int Zsz;                        /* input scalar */
    double *inMatrix;               /* 1xN input matrix */
    double *outMatrix;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first two arguments are scalar */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Number of clients must be a scalar.");
    }
    
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Number of demand shock types must be a scalar.");
    }
    
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","CDF vector must be type double.");
    }

    /* get the value of the scalar inputs  */
    client_number = mxGetScalar(prhs[0]);
    Zsz = mxGetScalar(prhs[1]);
    
    /* check that number of rows in second input argument is 1 */
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a column vector of with size number of demand shocks.");
    }
    
    /* create a pointer to the real data in the input matrix  */
    #if MX_HAS_INTERLEAVED_COMPLEX
    inMatrix = mxGetDoubles(prhs[2]);
    #else
    inMatrix = mxGetPr(prhs[2]);
    #endif

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,Zsz,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    outMatrix = mxGetDoubles(plhs[0]);
    #else
    outMatrix = mxGetPr(plhs[0]);
    #endif

    /* call the computational routine */
    getTypeNumbers(client_number,Zsz,outMatrix,inMatrix);
}