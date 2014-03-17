#include "GMM.h"
#include <math.h>

using namespace hiddenMarkovModel;

GMM::GMM(const double *data, const unsigned int dataSize, std::vector<GaussState> states):
	dataSize(dataSize),
	data(data),
	changeThreshold(1e-4),
	states(states)
{
	this->pi = array1D(this->states.size(), 1.0/(double) this->states.size());
	this->chi = array2D(this->dataSize, array1D(this->states.size(), 0));
}


GMM::~GMM(void)
{
}

void GMM::updateChi(){
	for (unsigned int state = 0; state < this->states.size(); state += 1){
		for (unsigned int t = 0; t < this->dataSize; t += 1){
			this->chi[t][state] = this->pi[state] * this->states[state].pdf(this->data[t]);
		}
	}
	for (unsigned int t = 0; t < this->dataSize; t += 1){
		normaliseArray(this->chi[t]);
	}
}
double GMM::updatePi(){
	double changeSum2 = 0;
	double sum2 = 0;
	
	for (unsigned int state = 0; state < this->states.size(); state += 1){
		double sum = 0;
		for (unsigned int t = 0; t < this->dataSize; t += 1){
			sum += this->chi[t][state];
		}
		double oldPi = this->pi[state];
		this->pi[state] = sum / this->dataSize;

		changeSum2 += (oldPi - this->pi[state]) * (oldPi - this->pi[state]);
		sum2 += this->pi[state] * this->pi[state];
	}

	return sqrt(changeSum2 / sum2);
}
void GMM::updateMu(){
	for (unsigned int state = 0; state < this->states.size(); state += 1){
		double sum = 0;
		for (unsigned int t = 0; t < this->dataSize; t += 1){
			sum += this->data[t] * this->chi[t][state];
		}
		this->states[state].mean = sum / this->pi[state] / this->dataSize;
	}
}
void GMM::updateSigma(){
	for (unsigned int state = 0; state < this->states.size(); state += 1){
		double sum = 0;
		for (unsigned int t = 0; t < this->dataSize; t += 1){
			double diff = (this->data[t] - this->states[state].mean);
			sum += diff * diff * this->chi[t][state];
		}
		this->states[state].std = sqrt(sum / this->pi[state] / this->dataSize);
	}
}

double GMM::iterate(){
	double changeRatio;
	this->updateChi();
	changeRatio = this->updatePi();
	this->updateMu();
	this->updateSigma();
	return changeRatio;
}

double GMM::run(){
	double changeRatio = this->changeThreshold + 1;
	while (changeRatio > this->changeThreshold){
		changeRatio = this->iterate();
	}
	return changeRatio;
}
