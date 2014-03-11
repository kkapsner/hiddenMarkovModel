#include <assert.h>
#include <math.h>
#include "Binner.h"

using namespace hiddenMarkovModel;

Binner::Binner(unsigned int size, double *dataToFit, unsigned long dataSize){
	assert(size > 0);
	assert(dataSize > 0);

	this->size = size;
	this->min = dataToFit[0];
	this->max = dataToFit[0];

	for (unsigned long i = 1; i < dataSize; i += 1){
		if (dataToFit[i] < min){
			this->min = dataToFit[i];
		}
		else if (dataToFit[i] > max){
			this->max = dataToFit[i];
		}
	}

	if (this->max == this->min){
		this->max += 1;
	}

	this->binSize = (this->max - this->min) / this->size;
};

unsigned int Binner::getBinIndex(double value) const{
	assert(value >= this->min);
	assert(value <= this->max);

	if (value == this->max){
		return this->size - 1;
	}
	else {
		return (int) floor((value - this->min) / this->binSize);
	}
};
unsigned int Binner::operator[](double value) const{
	return this->getBinIndex(value);
};

double Binner::getBinValue(unsigned int bin) const{
	assert(bin < this->size);

	return (0.5 + (double) bin) * this->binSize + this->min;
};
double Binner::operator()(unsigned int bin) const{
	return this->getBinValue(bin);
};

unsigned int Binner::getSize() const{
	return this->size;
};

double Binner::getMin() const{
	return this->min;
}

double Binner::getMax() const{
	return this->max;
}

double Binner::getBinSize() const{
	return this->binSize;
}