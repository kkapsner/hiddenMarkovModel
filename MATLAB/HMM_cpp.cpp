#include "mex.h"
#include <vector>
#include "../../../c++/hiddenMarkovModel/InitialEmissionProbability.h" //change this to the correct path!
#include "../../../c++/hiddenMarkovModel/GaussState.cpp" //change this to the correct path!
#include "../../../c++/hiddenMarkovModel/HMMConfiguration.cpp" //change this to the correct path!
#include "../../../c++/hiddenMarkovModel/HMM.cpp" //change this to the correct path!
#include "../../../c++/hiddenMarkovModel/Binner.cpp" //change this to the correct path!

#define DEFAULT_FALSE(name){\
    value = mxGetField(options, 0, #name);\
	configuration.name = (value != NULL && mxIsLogicalScalarTrue(value));\
}
#define DEFAULT_TRUE(name){\
    value = mxGetField(options, 0, #name);\
	configuration.name = (value == NULL || mxIsLogicalScalarTrue(value));\
}
#define DEFAULT_DOUBLE(name, defaultValue){\
    value = mxGetField(options, 0, #name);\
	if (value != NULL && mxIsDouble(value) && mxGetM(value) * mxGetN(value) == 1){\
		configuration.name = mxGetScalar(value);\
	}\
	else {\
		configuration.name = defaultValue;\
	}\
}
#define DEFAULT_INT(name, defaultValue){\
    value = mxGetField(options, 0, #name);\
	if (value != NULL && mxIsDouble(value) && mxGetM(value) * mxGetN(value) == 1){\
		configuration.name = (int) mxGetScalar(value);\
	}\
	else {\
		configuration.name = defaultValue;\
	}\
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *data, *transition, *gauss, *statesOut;
	unsigned int dataSize[2], statesCount, iterationCount;

	/* Check for proper number of arguments. */
	if(nrhs < 3) {
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:HMM:invalidNumInputs",
			"Three input arguments required."
		);
	} else if(nlhs > 4) {
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:HMM:maxlhs",
			"Too many output arguments."
		);
	}



	/* The input must be a noncomplex double.*/
	if(
		!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
		!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
	){
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:HMM:inputNotRealDouble",
			"Input must be a noncomplex double."
		);
	}

	dataSize[0] = (int) mxGetM(prhs[0]);
	dataSize[1] = (int) mxGetN(prhs[0]);
	statesCount = (int) mxGetM(prhs[1]);
	if (mxGetN(prhs[1]) != statesCount){
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:HMM:noSquareTransitionMatrix",
			"Transition matrix has to be MxM."
		);
	}

	if (mxGetM(prhs[2]) != statesCount){
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:HMM:notMatchingStatesCount",
			"Gauss definition matrix has a different states count than the transition matrix."
		);
	}
	if (mxGetN(prhs[2]) != 2){
		mexErrMsgIdAndTxt(
			"MATLAB:Fit:HMM:invalidGaussDefinition",
			"Gauss definition matrix has to be Mx2."
		);
	}

	/* Create matrix for the return argument. */
	plhs[0] = mxCreateDoubleMatrix(dataSize[0], dataSize[1], mxREAL);

	/* Assign pointers to each input and output. */
	data = mxGetPr(prhs[0]);
	transition = mxGetPr(prhs[1]);
	gauss = mxGetPr(prhs[2]);

	statesOut = mxGetPr(plhs[0]);

	using namespace hiddenMarkovModel;

	HMMConfiguration configuration;
	if (nrhs > 3){
		mxArray const *options = prhs[3];
		mxArray *value;
		
		value = mxGetField(options, 0, "verbose");
		if (value != NULL && mxIsLogicalScalarTrue(value)){
			configuration.verbose = true;
			
			DEFAULT_FALSE(verboseOutputEmission);
			DEFAULT_TRUE(verboseOutputTransition);
		}
		
		DEFAULT_DOUBLE(minSelfTransition, 0);
		DEFAULT_DOUBLE(minEmission, 1e-6);
		
		DEFAULT_TRUE(doEmissionUpdate);
		DEFAULT_TRUE(doTransitionUpdate);
		
		DEFAULT_INT(binningCount, 300);
		DEFAULT_INT(maxIterations, 100);
		DEFAULT_INT(abortStateChanges, 5);
	}

	std::vector<InitialEmissionProbability*> initStates(0);
	for (unsigned int i = 0; i < statesCount; i += 1){
		initStates.push_back(new GaussState(gauss[i], gauss[statesCount + i]));
	}

	HMM model (data, dataSize[0] * dataSize[1], initStates, configuration.binningCount);

	// delete state pointers
	initStates.clear();

	for (unsigned int i = 0; i < statesCount; i += 1){
		for (unsigned int j = 0; j < statesCount; j += 1){
			model.setTransition(transition[i + statesCount * j], i, j);
		}
	}
	model.autoSetSelfTransition();

	model.configuration = &configuration;

	model.run(iterationCount);

	std::vector<unsigned int> states (dataSize[0] * dataSize[1], 0);

	model.viterbi(states);

	for (unsigned int i = 0; i < dataSize[0] * dataSize[1]; i += 1){
		statesOut[i] = (double) states[i] + 1;
	}
	
	if (nlhs > 1){
		/* Create matrix for the transition output. */
		plhs[1] = mxCreateDoubleMatrix(statesCount, statesCount, mxREAL);
		double *transitionOut = mxGetPr(plhs[1]);
		for (unsigned int from = 0; from < statesCount; from += 1){
			for (unsigned int to = 0; to < statesCount; to += 1){
				transitionOut[from + to * statesCount] = model.getTransition(from, to);
			}
		}
		
		if (nlhs > 2){
			/* Create matrix for the emission output. */
			plhs[2] = mxCreateDoubleMatrix(statesCount, configuration.binningCount, mxREAL);
			double *emissionOut = mxGetPr(plhs[2]);
			for (unsigned int state = 0; state < statesCount; state += 1){
				for (unsigned int bin = 0; bin < configuration.binningCount; bin += 1){
					emissionOut[state + bin * statesCount] = model.getEmissionPropability(state, bin);
				}
			}
			
			if (nlhs > 3){
				/* Create matrix for the iteration count output. */
				plhs[3] = mxCreateDoubleScalar((double) iterationCount);
			}
		}
	}
}