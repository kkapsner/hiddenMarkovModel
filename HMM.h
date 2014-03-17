#pragma once
#include "typedefs.h"
#include "Binner.h"
#include "HMMConfiguration.h"
#include "InitialEmissionProbability.h"

// for verbose output
#include <iostream>

namespace hiddenMarkovModel{

	class HMM
	{
	private:
		
		// size of the data
		unsigned long dataSize;
		// the dataset to model
		double *data;
			
		// binner
		Binner *binner;
			
		// number of states
		unsigned int stateCount;
		
		// alphas of the forward algorithm
		array2D alpha;
		// scaling factor of the alphas
		array1D scaling;
		// betas of the backward algorithm
		array2D beta;
		// gamma[i][j] = probability to be in state j at time point i
		array2D gamma;
		// state[i] = most probable state at time point i
		std::vector<unsigned int> state;
		// xi[i][j][k] = probability to go from state j to state k at time point i
		array3D xi;
			
		// transition matrix. size stateCount x stateCount.
		// transition[i][j] = probability to go from state i to state j
		array2D transition;
		
		// emission probability. size stateCount x binner->getSize()
		// emission[i][binner->getBinIndex(value)] = probability that state i emits value
		array2D emission;
	public:
		// 
		const HMMConfiguration configuration;

		// constructor
		HMM(double *data, unsigned long dataSize, std::vector<InitialEmissionProbability*> states, HMMConfiguration configuration = HMMConfiguration());
		// destructor
		~HMM();
			
		// void setTransitionMatrix(double **transitions);
		void setTransition(double transition, unsigned int from, unsigned int to);
		void autoSetSelfTransition();
		void autoSetSelfTransition(unsigned int from);

		double getTransition(unsigned int from, unsigned int to);
			
		// void setEmissionProbability(double **emission);
		// void setEmissionPropability(double *emission, int state);
		void setEmissionPropability(double *mean, double *std);
		void setEmissionPropability(array1D mean, array1D std);
		void setEmissionPropability(double mean, double std, unsigned int state);

		double getEmissionPropability(unsigned int state, unsigned int bin);
			
		void updateAlpha();
		void updateBeta();
		/**
		 * updates the state probabilities for every time point and the most propable state
		 * 
		 * @return int number of changed states
		 */
		unsigned int updateGamma();
		void updateXi();
		void updateTransition();
		void updateEmission();
			
		unsigned int iterate();
		unsigned int run();
		unsigned int run(unsigned int &iterationCount);
			
		void viterbi(std::vector<unsigned int> &states);

		// verbose output functions
		void outputEmission(std::ostream &out);
		void outputTransition(std::ostream &out);
		void outputState(std::ostream &out);
	};
};

