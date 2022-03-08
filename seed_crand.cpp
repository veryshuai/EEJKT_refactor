#include "mex.h"
#include <cstdlib>
#include <ctime>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{  
    
    int rand_seed;

    rand_seed = mxGetScalar(prhs[0]);
    srand(rand_seed);
    return;
    
}
