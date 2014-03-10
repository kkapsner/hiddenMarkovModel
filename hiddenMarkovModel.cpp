// hiddenMarkovModel.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//

#include "stdafx.h"
#include <iostream>


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

	HMMConfiguration configuration;
	
	std::cout << "create model" << std::endl;
	HMM model (data, dataSize, states, 100);

	std::cout << "set transition" << std::endl;
	for (int i = 0; i < 2; i += 1){
		for (int j = 0; j < 2; j += 1){
			model.setTransition(0.5, i, j);
		}
	}
  
	model.configuration = &configuration;
	configuration.verbose = true;
	configuration.pauseAfterIteration = true;

	model.run();

	std::vector<unsigned int> fittedStates (dataSize, 0);

	model.viterbi(fittedStates);

	delete data;
	states.clear();

	/*
	array1D test (20, 0);
	double c = 0;
	for (array1D::iterator iter = test.begin(); iter < test.end(); iter += 1){
		*iter = c;
		c += 2;
	}
	rescaleArray(test, 2);
	for (array1D::iterator iter = test.begin(); iter < test.end(); iter += 1){
		std::cout << *iter << std::endl;
	}

	std::cout << std::endl << normaliseArray(test) << std::endl << std::endl;
	
	for (array1D::iterator iter = test.begin(); iter < test.end(); iter += 1){
		std::cout << *iter << std::endl;
	}*/

	std::cout << "Hit enter to exit." << std::endl;
	char in[1];
	std::cin.read(in, 1);
	return 0;
}

