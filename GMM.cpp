#include "GMM.h"
#include <math.h>
#include <algorithm>
#include <iostream>

using namespace hiddenMarkovModel;

GMM::GMM(array1D data, std::vector<GaussState> states):
	data(data),
	changeThreshold(1e-4),
	states(states)
{
	this->pi = array1D(this->states.size(), 1.0/(double) this->states.size());
	this->chi = array2D(this->data.size(), array1D(this->states.size(), 0));
}

GMM::GMM(array1D data, std::vector<GaussState> states, array1D pi):
	data(data),
	changeThreshold(1e-4),
	states(states),
	pi(pi)
{
	this->chi = array2D(this->data.size(), array1D(this->states.size(), 0));
}


GMM::~GMM(void)
{
}

void GMM::initialiseStates(){
	unsigned int dataSize = this->data.size();
	unsigned int statesCount = this->states.size();
	// sort states;
	std::sort(this->data.begin(), this->data.end());

	std::vector<unsigned int> count(statesCount, 0);
	array1D mu(statesCount, 0);
	array1D sigma(statesCount, 0);

	for (unsigned int t = 0; t < dataSize; t += 1){
		unsigned int state = (t * statesCount) / dataSize;
		count[state] += 1;
		mu[state] += this->data[t];
	}

	for (unsigned int state = 0; state < statesCount; state += 1){
		this->pi[state] = (double) count[state] / (double) dataSize;
		this->states[state].mean = mu[state] /(double) count[state];
	}
	
	for (unsigned int t = 0; t < dataSize; t += 1){
		unsigned int state = (t * statesCount) / dataSize;
		double diff = this->data[t] - this->states[state].mean;
		sigma[state] += diff * diff;
	}

	for (unsigned int state = 0; state < statesCount; state += 1){
		this->states[state].std = sqrt(sigma[state] / (double) count[state]);
	}
}

void GMM::updateChi(){
	unsigned int dataSize = this->data.size();
	unsigned int statesCount = this->states.size();
	for (unsigned int t = 0; t < dataSize; t += 1){
		for (unsigned int state = 0; state < statesCount; state += 1){
			this->chi[t][state] = this->pi[state] * this->states[state].pdf(this->data[t]);
		}
		normaliseArray(this->chi[t]);
	}
}
double GMM::updatePi(){
	unsigned int dataSize = this->data.size();
	unsigned int statesCount = this->states.size();
	double changeSum2 = 0;
	double sum2 = 0;
	
	for (unsigned int state = 0; state < statesCount; state += 1){
		double sum = 0;
		for (unsigned int t = 0; t < dataSize; t += 1){
			sum += this->chi[t][state];
		}
		double oldPi = this->pi[state];
		this->pi[state] = sum / dataSize;

		changeSum2 += (oldPi - this->pi[state]) * (oldPi - this->pi[state]);
		sum2 += this->pi[state] * this->pi[state];
	}

	return sqrt(changeSum2 / sum2);
}
double GMM::updateMu(){
	unsigned int dataSize = this->data.size();
	unsigned int statesCount = this->states.size();
	double changeSum2 = 0;
	double sum2 = 0;

	for (unsigned int state = 0; state < statesCount; state += 1){
		double sum = 0;
		for (unsigned int t = 0; t < dataSize; t += 1){
			sum += this->data[t] * this->chi[t][state];
		}

		double oldMean = this->states[state].mean;
		this->states[state].mean = sum / this->pi[state] / dataSize;

		changeSum2 += (oldMean - this->states[state].mean) * (oldMean - this->states[state].mean);
		sum2 += this->states[state].mean * this->states[state].mean;
	}

	return sqrt(changeSum2 / sum2);
}
double GMM::updateSigma(){
	unsigned int dataSize = this->data.size();
	unsigned int statesCount = this->states.size();
	double changeSum2 = 0;
	double sum2 = 0;

	for (unsigned int state = 0; state < statesCount; state += 1){
		double sum = 0;
		for (unsigned int t = 0; t < dataSize; t += 1){
			double diff = (this->data[t] - this->states[state].mean);
			sum += diff * diff * this->chi[t][state];
		}
		double oldStd = this->states[state].std;
		this->states[state].std = sqrt(sum / this->pi[state] / dataSize);

		changeSum2 += (oldStd - this->states[state].std) * (oldStd - this->states[state].std);
		sum2 += this->states[state].std * this->states[state].std;
	}

	return sqrt(changeSum2 / sum2);
}

double GMM::iterate(){
	double changeRatio = 0;
	this->updateChi();
	changeRatio += this->updatePi();
	changeRatio += this->updateMu();
	changeRatio += this->updateSigma();
	return changeRatio / 3;
}

double GMM::run(){
	unsigned int i;
	return this->run(i);
}

double GMM::run(unsigned int &iterationCount){
	iterationCount = 0;
	double changeRatio = this->changeThreshold + 1;
	while (changeRatio > this->changeThreshold){
		iterationCount += 1;
		std::cout << "Iteration " << iterationCount;
		changeRatio = this->iterate();
		std::cout << " (" << changeRatio << ")" << std::endl;
		
		for (unsigned int state = 0; state < this->states.size(); state += 1){
			std::cout << "State " << (state + 1) << ": " <<
				"pi: " << this->pi[state] << ", "
				"mu: " << this->states[state].mean << ", "
				"sigma: " << this->states[state].std << std::endl;
		}
		std::cout << std::endl;
	}
	return changeRatio;
}
