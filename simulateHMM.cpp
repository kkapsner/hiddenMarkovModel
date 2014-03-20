#pragma once
#include "typedefs.h"
#include "GaussState.h"
#include <random>

namespace hiddenMarkovModel{
	void simulateHMM(
		array1D &data,
		std::vector<unsigned int> &state,
		const array2D transition,
		const std::vector<GaussState> states,
		const unsigned int dataSize
	){
		const unsigned int statesCount = states.size();

		data.resize(dataSize);
		state.resize(dataSize);
		state[0] = std::rand() % statesCount;
		data[0] = states[state[0]].random();

		for (unsigned int t = 1; t < dataSize; t += 1){
			double r = (double) rand() / (double) RAND_MAX;
			unsigned int s;
			for (s = 0; s < statesCount; s += 1){
				r -= transition[state[t - 1]][s];
				if (r <= 0){
					break;
				}
			}
			if (s == 3){
				s = 2;
			}
			state[t] = s;
			data[t] = states[state[t]].random();
		}

	};
}