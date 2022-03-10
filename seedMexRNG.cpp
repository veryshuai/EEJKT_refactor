/*
 * seedMexRNG.cpp 
 * 
 * this function seeds the random number generator for mex files called subsequently
 *
 *      call like this:
 *      seedMexRNG(SEED INTEGER);
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include <random> 
#include <iostream> 
#include <stdlib.h>
#include "math.h" 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//this is the link function
{

    double rand; //declare random number

    //seed the random number generator 
    srand(mxGetScalar(prhs[0]));

    //print a rand to test
    //rand = std::rand();
    // mexPrintf("A random number: %4.0f \n",rand);

}
