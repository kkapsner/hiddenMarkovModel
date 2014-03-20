#include "GaussState.h"
#include <math.h>
#include <random>

using namespace hiddenMarkovModel;

GaussState::GaussState(void):
	mean(0),
	std(1)
	{}

GaussState::GaussState(double mean, double std):
	mean(mean),
	std(std)
	{}


GaussState::~GaussState(void){
}

double GaussState::pdf(double x) const{
    static const double inv_sqrt_2pi = 0.3989422804014327;
	double a = (x - this->mean) / this->std;
	return inv_sqrt_2pi / this->std * exp(- a * a / 2.0);
}

double GaussState::random() const{
	//return this->mean;
	std::random_device rd;
    std::mt19937 gen(rd());
	std::normal_distribution<> d(this->mean, this->std);
	return d(gen);/**/
}
