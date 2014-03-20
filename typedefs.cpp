#pragma once
#include "typedefs.h"

using namespace hiddenMarkovModel;


void hiddenMarkovModel::rescaleArray(array1D &arr, double factor){
	for (unsigned int i = 0, end = arr.size(); i < end; i += 1){
		arr[i] /= factor;
	}
};

double hiddenMarkovModel::normaliseArray(array1D &arr){
	double sum = 0;
	for (unsigned int i = 0, end = arr.size(); i < end; i += 1){
		sum += arr[i];
	}
	rescaleArray(arr, sum);
	return sum;
};