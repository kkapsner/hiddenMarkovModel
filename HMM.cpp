#include "HMM.h"
#include <assert.h>
#include "Binner.h"
#include "InitialEmissionProbability.h"
#include <iostream>

using namespace hiddenMarkovModel;

double hiddenMarkovModel::normaliseArray(array1D &arr){
	double sum = 0;
	for (array1D::iterator iter = arr.begin(); iter < arr.end(); iter++){
		sum += *iter;
	}
	rescaleArray(arr, sum);
	return sum;
}
void hiddenMarkovModel::rescaleArray(array1D &arr, double factor){
	for (array1D::iterator iter = arr.begin(); iter < arr.end(); iter++){
		*iter /= factor;
	}
}

HMM::HMM(double *data, unsigned long dataSize, std::vector<InitialEmissionProbability*> states, unsigned int binnerSize){
	unsigned int stateCount = states.size();
	assert(dataSize > 0);
	assert(stateCount > 0);
	assert(binnerSize > 0);

	this->data = data;
	this->dataSize = dataSize;
	this->stateCount = stateCount;
	this->binner = new Binner(binnerSize, data, dataSize);

	// initialise 1D-arrays
	this->scaling = array1D(this->dataSize, 1);
	this->state = std::vector<unsigned int>(this->dataSize, stateCount);

	// initialise 2D-arrays
	
	this->alpha = array2D(this->dataSize, array1D(this->stateCount, 1.0 / (double) this->stateCount));
	this->beta = array2D(this->dataSize, array1D(this->stateCount, 1.0 / (double) this->stateCount));
	this->gamma = array2D(this->dataSize, array1D(this->stateCount, 1.0 / (double) this->stateCount));

	// fill emission probability according to the given initial 
	this->emission = array2D(this->stateCount, array1D(binnerSize, 0));
	for (unsigned int bin = 0; bin < this->binner->getSize(); bin += 1){
		double x = this->binner->getBinValue(bin);
		for (unsigned int i = 0; i < this->stateCount; i += 1){
			this->emission[i][bin] = states[i]->pdf(x);
		}
	}
	for (unsigned int i = 0; i < this->stateCount; i += 1){
		normaliseArray(this->emission[i]);
	}

	this->transition = array2D(this->stateCount, array1D(this->stateCount, 0));

	this->xi = array3D(this->dataSize, array2D(this->stateCount, array1D(this->stateCount, 1.0 / (double) this->stateCount)));
}

HMM::~HMM(){
	delete this->binner;
}


void HMM::setEmissionPropability(double *mean, double *std){
	for (unsigned int i = 0; i < this->stateCount; i += 1){
		this->setEmissionPropability(mean[i], std[i], i);
	}
}
void HMM::setEmissionPropability(array1D mean, array1D std){
	assert(mean.size() >= this->stateCount);
	assert(std.size() >= this->stateCount);
	for (unsigned int i = 0; i < this->stateCount; i += 1){
		this->setEmissionPropability(mean[i], std[i], i);
	}
}
void HMM::setEmissionPropability(double mean, double std, unsigned int state){
	// translate mean and std in binned space
	mean = (mean - this->binner->getMin()) / (double) this->binner->getBinSize();
	std = std / (double) this->binner->getBinSize();

	for (unsigned int i = 0; i < this->binner->getSize(); i += 1){
		double x = ((double) i - mean) / std;
		this->emission[state][i] = exp(- x * x);
	}
	normaliseArray(this->emission[state]);
}

double HMM::getEmissionPropability(unsigned int state, unsigned int bin){
	return this->emission[state][bin];
}

void HMM::setTransition(double transition, unsigned int from, unsigned int to){
	this->transition[from][to] = transition;
}

void HMM::autoSetSelfTransition(){
	for (unsigned int from = 0; from < this->stateCount; from += 1){
		this->autoSetSelfTransition(from);
	}
}
void HMM::autoSetSelfTransition(unsigned int from){
	double sum = 0;
	for (unsigned int to = 0; to < this->stateCount; to += 1){
		if (to != from){
			sum += this->transition[from][to];
		}
	}
	assert(sum <= 1);
	this->transition[from][from] = 1 - sum;
}

double HMM::getTransition(unsigned int from, unsigned int to){
	return this->transition[from][to];
}

void HMM::updateAlpha(){
	// first time point
	for (unsigned int i = 0; i < this->stateCount; i += 1){
		this->alpha[0][i] = this->gamma[0][i] * this->emission[i][(*this->binner)[this->data[0]]];
	}
	this->scaling[0] = normaliseArray(this->alpha[0]);

	for (unsigned int i = 1; i < this->dataSize; i += 1){
		for (unsigned int j = 0; j < this->stateCount; j += 1){
			this->alpha[i][j] = 0;
			for (unsigned int k = 0; k < this->stateCount; k += 1){
				this->alpha[i][j] += this->alpha[i-1][k] * this->transition[k][j] * this->emission[j][(*this->binner)[this->data[i]]];
			}
		}
		this->scaling[i] = normaliseArray(this->alpha[i]);
	}
}

void HMM::updateBeta(){
	for (unsigned int i = 0; i < this->stateCount; i += 1){
		this->beta[this->dataSize - 1][i] = 1;
	}
	for (int i = this->dataSize - 2; i >= 0; i -= 1){
		for (unsigned int j = 0; j < this->stateCount; j += 1){
			this->beta[i][j] = 0;
			for (unsigned int k = 0; k < this->stateCount; k += 1){
				this->beta[i][j] += this->beta[i+1][k] * this->transition[j][k] * this->emission[k][(*this->binner)[this->data[i+1]]];
			}
			
		}
		rescaleArray(this->beta[i], this->scaling[i]);
	}
}

unsigned int HMM::updateGamma(){
	unsigned int changeCount = 0;
	for (unsigned int i=0; i < this->dataSize; i += 1){
		double denominator = 0;
		for (unsigned int j = 0; j < this->stateCount; j += 1){
			denominator += this->alpha[i][j] * this->beta[i][j];
		}
		double maxGamma = -1;
		unsigned int maxState = this->stateCount;
		for (unsigned int j = 0; j < this->stateCount; j += 1){
			this->gamma[i][j] = this->alpha[i][j] * this->beta[i][j] / denominator;
			if (this->gamma[i][j] > maxGamma){
				maxState = j;
				maxGamma = this->gamma[i][j];
			}
		}

		if (maxState != this->state[i]){
			changeCount += 1;
			this->state[i] = maxState;
		}
	}
	return changeCount;
}

void HMM::updateXi(){
	for (unsigned int t = 0; t < this->dataSize - 1; t += 1){
		for (unsigned int i = 0; i < this->stateCount; i += 1){
			for (unsigned int j = 0; j < this->stateCount; j += 1){
				this->xi[t][i][j] = (
					this->gamma[t][i] *
					this->transition[i][j] *
					this->emission[j][(*this->binner)[this->data[t+1]]] *
					this->scaling[t+1] *
					this->beta[t+1][j] /
					this->beta[t][i]
				);
			}
		}
	}
}

void HMM::updateEmission(){
	for (unsigned int i = 0; i < this->stateCount; i += 1){
		for (unsigned int j = 0; j < this->binner->getSize(); j += 1){
			this->emission[i][j] = 0;
		}
	}

	for (unsigned int t = 0; t < this->dataSize; t += 1){
		unsigned int currIndex = (*this->binner)[this->data[t]];

		for (unsigned int i = 0; i < this->stateCount; i += 1){
			this->emission[i][currIndex] += this->gamma[t][i];
		}
	}
	for (unsigned int i = 0; i < this->stateCount; i += 1){
		normaliseArray(this->emission[i]);

		bool renorm = false;
		for (unsigned int j = 0; j < this->binner->getSize(); j += 1){
			if (this->emission[i][j] < this->configuration->minEmission){
				this->emission[i][j] = this->configuration->minEmission;
				renorm = true;
			}
		}
		if (renorm){
			normaliseArray(this->emission[i]);
		}
	}
}

void HMM::updateTransition(){

	for (unsigned int i = 0; i < this->stateCount; i += 1){
		for (unsigned int j = 0; j < this->stateCount; j += 1){
			this->transition[i][j] = 0;
			for (unsigned int t = 0; t < this->dataSize - 1; t += 1){
				this->transition[i][j] += this->xi[t][i][j];
			}
		}
		normaliseArray(this->transition[i]);

		// restrict transition rates

		if (this->transition[i][i] < this->configuration->minSelfTransition){
			this->transition[i][i] = this->configuration->minSelfTransition * (1 + this->transition[i][i]) / (1 - this->configuration->minSelfTransition);
			normaliseArray(this->transition[i]);
		}
	}
}

unsigned int HMM::iterate(){
	int changeCount;
	this->updateAlpha();
	this->updateBeta();
	changeCount = this->updateGamma();
	this->updateXi();
	if (this->configuration->doEmissionUpdate){
		this->updateEmission();
	}
	if (this->configuration->doTransitionUpdate){
		this->updateTransition();
	}

	return changeCount;
}

unsigned int HMM::run(){
	unsigned int i;
	return this->run(i);
}

unsigned int HMM::run(unsigned int &i){
	unsigned int changeCount = 0;
	for (i = 0; i < this->configuration->maxIterations; i += 1){
		if (this->configuration->verbose){
			std::cout << std::endl << "start iteration " << (i + 1) << ":" << std::endl; 
		}
		changeCount = this->iterate();
		if (this->configuration->verbose){
			std::cout << "changed states: " << changeCount << std::endl;
			if (this->configuration->verboseOutputEmission){
				std::cout << "emission probabilities:" << std::endl;
				this->outputEmission(std::cout);
			}
			if (this->configuration->verboseOutputTransition){
				std::cout << "transition probabilities:" << std::endl;
				this->outputTransition(std::cout);
			}
		}
		if (changeCount <= this->configuration->abortStateChanges){
			i += 1;
			break;
		}

		if (this->configuration->pauseAfterIteration){
			std::cout << "PAUSE" << std::endl << std::endl << "Hit enter to resume." << std::endl;
			char in[1];
			std::cin.read(in, 1);
		}
	}
	return changeCount;
}

void HMM::viterbi(std::vector<unsigned int> &states){
	states.resize(this->dataSize);
	/*for (unsigned int i = 0; i < this->dataSize; i += 1){
		states[i] = this->state[i];
	}
	*/

	array2D logEmission (this->stateCount, array1D (this->binner->getSize(), 0));
	for (unsigned int state = 0; state < this->stateCount; state += 1){
		for (unsigned int bin = 0; bin < this->binner->getSize(); bin += 1){
			logEmission[state][bin] = log(this->emission[state][bin]);
		}
	}
	array2D logTransition (this->stateCount, array1D (this->stateCount, 0));
	for (unsigned int from = 0; from < this->stateCount; from += 1){
		for (unsigned int to = 0; to < this->stateCount; to += 1){
			logTransition[from][to] = log(this->transition[from][to]);
		}
	}
	
	array2D v (this->dataSize, array1D (this->stateCount, 0));
	std::vector<std::vector<unsigned int>> psi (this->dataSize, std::vector<unsigned int> (this->stateCount, 0));
	for (unsigned int state = 0; state < this->stateCount; state += 1){
		v[0][state] = log(this->gamma[0][state]) + logEmission[state][(*this->binner)[this->data[0]]];
		psi[0][state] = 0;
	}

	for (unsigned int t = 1; t < this->dataSize; t += 1){
		for (unsigned int state = 0; state < this->stateCount; state += 1){
			double maxValue = v[t-1][0] + logTransition[0][state] + logEmission[state][(*this->binner)[this->data[t]]];

			for (unsigned int lastState = 1; lastState < this->stateCount; lastState += 1){
				double value = v[t-1][lastState] + logTransition[lastState][state] + logEmission[state][(*this->binner)[this->data[t]]];
				if (value > maxValue){
					maxValue = value;
				}
			}
			v[t][state] = maxValue;
		}

		for (unsigned int state = 0; state < this->stateCount; state += 1){
			unsigned int maxIndex = 0;
			double maxValue = v[t][0] + logTransition[0][state];
			for (unsigned int from = 1; from < this->stateCount; from += 1){
				double value = v[t][from] + logTransition[from][state];
				if (value > maxValue){
					maxValue = value;
					maxIndex = from;
				}
			}
			psi[t][state] = maxIndex;
		}
	}

	unsigned int maxEndState = 0;
	double value = v[this->dataSize - 1][0];
	for (unsigned int state = 1; state < this->stateCount; state += 1){
		if (v[this->dataSize - 1][state] > value){
			value = v[this->dataSize - 1][state];
			maxEndState = state;
		}
	}

	states[this->dataSize - 1] = maxEndState;
	for (int t = this->dataSize - 2; t >= 0; t -= 1){
		states[t] = psi[t + 1][states[t + 1]];
	}
}

void HMM::outputEmission(std::ostream &out){
	for (unsigned int bin = 0; bin < this->binner->getSize(); bin += 1){
		out << bin << ":\t";
		for (unsigned int state = 0; state < this->stateCount; state += 1){
			out << this->emission[state][bin] << "\t";
		}
		out << std::endl;
	}
}
void HMM::outputTransition(std::ostream &out){
	for (unsigned int from = 0; from < this->stateCount; from += 1){
		for (unsigned int to = 0; to < this->stateCount; to += 1){
			out << this->transition[from][to] << "\t";
		}
		out << std::endl;
	}
}
void HMM::outputState(std::ostream &out){
	out << this->state[0];
	for (unsigned int t = 0; t < this->dataSize; t += 1){
		out << ", " << this->state[t];
	}
	out << std::endl;
}