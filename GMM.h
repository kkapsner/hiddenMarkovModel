#pragma once
#include "typedefs.h"
#include "GaussState.h"

namespace hiddenMarkovModel{
	class GMM{
		private:
			array1D data;
			array1D pi;
			array2D chi;
		public:
			std::vector<GaussState> states;
			double changeThreshold;
			
			GMM(array1D data, std::vector<GaussState> states);
			GMM(array1D data, std::vector<GaussState> states, array1D pi);
			~GMM(void);

			void initialiseStates();

			void updateChi();
			double updatePi();
			double updateMu();
			double updateSigma();

			double iterate();
			double run();
			double run(unsigned int &iterationCount);
	};
}

