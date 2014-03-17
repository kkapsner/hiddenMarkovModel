#pragma once
#include "typedefs.h"
#include "GaussState.h"

namespace hiddenMarkovModel{
	class GMM{
		private:
			const double *data;
			const unsigned int dataSize;
			array1D pi;
			array2D chi;
		public:
			std::vector<GaussState> states;
			double changeThreshold;

			GMM(const double *data, const unsigned int dataSize, std::vector<GaussState> states);
			~GMM(void);

			void updateChi();
			double updatePi();
			void updateMu();
			void updateSigma();

			double iterate();
			double run();
	};
}

