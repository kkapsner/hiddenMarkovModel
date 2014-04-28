// hiddenMarkovModel.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>


int _tmain(int argc, _TCHAR* argv[])
{
	// only tests up to now
	using namespace hiddenMarkovModel;

	int dataSize = 5000;

	double *data = new double[dataSize];
	for (int i = 0; i < dataSize; i += 1){
		data[i] = 1.0/(double)(i + 1);
	}
	
	std::vector<InitialEmissionProbability*> states(0);
	states.push_back(new GaussState(0, 1));
	states.push_back(new GaussState(1, 1));

	// HMMConfiguration configuration;
	// configuration.binningCount = 100;
	// configuration.verbose = true;
	// configuration.pauseAfterIteration = true;
	HMMConfiguration configuration = HMMConfiguration::fromFile(std::ifstream("test.json"));
	
	std::cout << "create model" << std::endl;
	HMM model (data, dataSize, states, configuration);

	array1D binnerRange(2, 0);
	model.getBinningRange(binnerRange);
	std::cout << "binning range: " << binnerRange[0] << " to " << binnerRange[1]<< std::endl;

	std::cout << "set transition" << std::endl;
	for (int i = 0; i < 2; i += 1){
		for (int j = 0; j < 2; j += 1){
			model.setTransition(0.5, i, j);
		}
	}

	model.run();

	std::vector<unsigned int> fittedStates (dataSize, 0);

	model.viterbi(fittedStates);

	delete data;
	states.clear();


	std::cout << "Hit enter to exit." << std::endl;
	char in[1];
	std::cin.read(in, 1);
	return 0;
}

// hiddenMarkovModel.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//
/*
#include "stdafx.h"
#include "GMM.h"
#include <iostream>
#include <fstream>


int _tmain(int argc, _TCHAR* argv[])
{
	// only tests up to now
	using namespace hiddenMarkovModel;

	unsigned int dataSize = 100000;
	array1D data(dataSize);
	std::vector<unsigned int> state(dataSize);

	array2D transition(3, array1D(3, 0));
	transition[0][0] = 0.8;
	transition[0][1] = 0.1;
	transition[0][2] = 0.1;
	transition[1][0] = 0.2;
	transition[1][1] = 0.6;
	transition[1][2] = 0.2;
	transition[2][0] = 0.1;
	transition[2][1] = 0.1;
	transition[2][2] = 0.8;
	
	std::vector<GaussState> states(0);
	states.push_back(GaussState(0, 20));
	states.push_back(GaussState(50, 20));
	states.push_back(GaussState(100, 20));
	
	std::cout << "simulate data" << std::endl;
	simulateHMM(data, state, transition, states, dataSize);
	
	std::cout << "create model" << std::endl;
	GMM model (data, states);

	model.run();

	std::cout << "Hit enter to exit." << std::endl;
	char in[1];
	std::cin.read(in, 1);
	return 0;
}
*/
