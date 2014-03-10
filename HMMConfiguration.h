#pragma once

namespace hiddenMarkovModel {
	class HMMConfiguration{
	public:
		bool verbose;
		bool verboseOutputTransition;
		bool verboseOutputEmission;
		bool pauseAfterIteration;

		double minSelfTransition;
		double minEmission;

		bool doEmissionUpdate;
		bool doTransitionUpdate;
		unsigned int binningCount;
		unsigned int maxIterations;
		unsigned int abortStateChanges;
		HMMConfiguration();
	};
};

