#pragma once

#include "InitialEmissionProbability.h"

namespace hiddenMarkovModel{
	class GaussState:
		public InitialEmissionProbability
	{
	public:
		double mean;
		double std;
		GaussState();
		GaussState(double mean, double std);
		~GaussState();

		virtual double pdf(double x);
		virtual double random();
	};
};

