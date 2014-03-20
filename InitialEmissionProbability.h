#pragma once
namespace hiddenMarkovModel{
	class InitialEmissionProbability{
		public:
			virtual double pdf(double x) const = 0;
			virtual double random() const = 0;
	};
}
