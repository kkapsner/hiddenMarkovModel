#include "HMMConfiguration.h"

using namespace hiddenMarkovModel;

HMMConfiguration::HMMConfiguration():
	verbose(false),
	verboseOutputTransition(true),
	verboseOutputEmission(false),
	
	pauseAfterIteration(false),

	minSelfTransition(0),
	minEmission(1e-6),

	doEmissionUpdate(true),
	doTransitionUpdate(true),
	binningCount(300),
	maxIterations(100),
	abortStateChanges(5)
	{}
