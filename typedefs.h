#pragma once
#include <vector>

namespace hiddenMarkovModel{
	typedef std::vector<double> array1D;
	typedef std::vector<array1D> array2D;
	typedef std::vector<array2D> array3D;
	

	void rescaleArray(array1D &arr, double factor);

	double normaliseArray(array1D &arr);
};