#pragma once
#include <vector>
#include "GaussState.h"

namespace hiddenMarkovModel{
	typedef std::vector<double> array1D;
	typedef std::vector<array1D> array2D;
	typedef std::vector<array2D> array3D;
	

	void rescaleArray(array1D &arr, double factor);

	double normaliseArray(array1D &arr);

	void simulateHMM(
		array1D &data,
		std::vector<unsigned int> &state,
		const array2D transition,
		const std::vector<GaussState> states,
		const unsigned int dataSize
	);
};