#pragma once

namespace hiddenMarkovModel{
	class Binner
	{
		private:
			// number of bins
			unsigned int size;
			// lower boundary
			double min;
			// upper boundary
			double max;
			// bin size
			double binSize;
			
		public:
			// constructor
			Binner(unsigned int size, double* dataToFit, unsigned long dataSize);

			/**
			 * returns the bin to a given value. Throws an exception if the value is out of bound.
			 *
			 * @param double value
			 * @return unsigned int
			 */
			unsigned int getBinIndex(double value) const;
			unsigned int operator [](double value) const;
			
			/**
			 * returns the middle value of the bin. Throws an exception if the bin is out of bound.
			 *
			 * @param unsigned int bin
			 * @return double
			 */
			double getBinValue(unsigned int bin) const;
			double operator ()(unsigned int bin) const;
			
			/**
			 * returns the size of the binner
			 *
			 * @return int
			 */
			unsigned int getSize() const;
			double getMin() const;
			double getMax() const;
			double getBinSize() const;
	};
}

