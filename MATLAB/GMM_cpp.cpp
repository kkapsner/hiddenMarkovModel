#include "mex.h"
#include "../typedefs.cpp"
#include <vector>
#include "../GaussState.cpp"
#include "../GMM.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *data, *mu, *sigma;
	unsigned int dataSize[2], dataPointCount, dataStart, dataEnd, statesCount, iterationCount;

	/* Check for proper number of arguments. */
	if(nrhs < 2) {
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:GMM:invalidNumInputs",
			"Two input arguments required."
		);
	} else if(nlhs < 2) {
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:GMM:maxlhs",
			"Too few output arguments."
		);
	} else if(nlhs > 4) {
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:GMM:maxlhs",
			"Too many output arguments."
		);
	}



	/* The input must be a noncomplex double.*/
	if(
		!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
		!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
	){
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:GMM:inputNotRealDouble",
			"Input must be a noncomplex double."
		);
	}

	dataSize[0] = (int) mxGetM(prhs[0]);
	dataSize[1] = (int) mxGetN(prhs[0]);
	dataPointCount = dataSize[0] * dataSize[1];
	if (mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1){
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:GMM:noScalarStatesCount",
			"States count has to be scalar."
		);
	}
	statesCount = mxGetScalar(prhs[1]);

	/* Create matrix for the return argument. */
	plhs[0] = mxCreateDoubleMatrix(statesCount, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(statesCount, 1, mxREAL);

	/* Assign pointers to each input and output. */
	data = mxGetPr(prhs[0]);

	mu = mxGetPr(plhs[0]);
	sigma = mxGetPr(plhs[1]);

	using namespace hiddenMarkovModel;

	std::vector<GaussState> states(statesCount);
	
	// check for NaN at the start
	dataStart = 0;
	for (unsigned int i = 0; i < dataPointCount; i += 1){
		if (!mxIsNaN(data[i])){
			dataStart = i;
			break;
		}
	}
	// check for NaN at the end
	for (dataEnd = dataStart; dataEnd < dataPointCount; dataEnd += 1){
		if (mxIsNaN(data[dataEnd])){
			break;
		}
	}

	GMM model (std::vector<double> (data + dataStart, data + dataEnd), states);

	model.initialiseStates();

	model.run(iterationCount);

	for (unsigned int state = 0; state < statesCount; state += 1){
		mu[state] = model.states[state].mean;
		sigma[state] = model.states[state].std;
	}
	
	if (nlhs > 2){
		/* Create matrix for the gaussian mixture output. */
		plhs[2] = mxCreateDoubleMatrix(statesCount, 1, mxREAL);
		double *mixture;
		mixture = mxGetPr(plhs[2]);
		for (unsigned int state = 0; state < statesCount; state += 1){
			model.getPi(mixture);
		}
		
	}
	if (nlhs > 3){
		/* Create matrix for the iteration count output. */
		plhs[3] = mxCreateDoubleScalar((double) iterationCount);
	}
}