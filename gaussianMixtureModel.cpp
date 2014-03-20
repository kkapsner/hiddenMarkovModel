// hiddenMarkovModel.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//

#include "stdafx.h"
#include "simulateHMM.cpp"
#include "GMM.h"
#include <iostream>
#include <fstream>


int _tmain(int argc, _TCHAR* argv[])
{
	// only tests up to now
	using namespace hiddenMarkovModel;

	array1D data;
	std::vector<unsigned int> state;

	array2D transition(3, array1D(3, 0));
	transition[0][0] = 0.8;
	transition[0][1] = 0.1;
	transition[0][2] = 0.1;
	transition[1][0] = 0.6;
	transition[1][1] = 0.2;
	transition[1][2] = 0.2;
	transition[2][0] = 0.8;
	transition[2][1] = 0.1;
	transition[2][2] = 0.1;
	
	std::vector<GaussState> states(0);
	states.push_back(GaussState(0, 20));
	states.push_back(GaussState(50, 20));
	states.push_back(GaussState(100, 20));
	int dataSize = 1000000;

	simulateHMM(data, state, transition, states, dataSize);
	
	std::cout << "create model" << std::endl;
	GMM model (data, states);

	model.run();

	std::cout << "Hit enter to exit." << std::endl;
	char in[1];
	std::cin.read(in, 1);
	return 0;
}

