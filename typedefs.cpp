#pragma once
#include "typedefs.h"

using namespace hiddenMarkovModel;


void hiddenMarkovModel::rescaleArray(array1D &arr, double factor){
	for (array1D::iterator iter = arr.begin(); iter < arr.end(); iter++){
		*iter /= factor;
	}
};

double hiddenMarkovModel::normaliseArray(array1D &arr){
	double sum = 0;
	for (array1D::iterator iter = arr.begin(); iter < arr.end(); iter++){
		sum += *iter;
	}
	rescaleArray(arr, sum);
	return sum;
};