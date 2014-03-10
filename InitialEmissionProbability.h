#pragma once
namespace hiddenMarkovModel{
	class InitialEmissionProbability{
		public:
			virtual double pdf(double x) = 0;
			virtual double random() = 0;
	};
}
