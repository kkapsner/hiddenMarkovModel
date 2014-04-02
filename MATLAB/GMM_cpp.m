%GMM_CPP performs a gaussian mixture model (GMM) analysis on the given dataset
%
%Parameter:
%   data: vector of the data to be modeled with a GMM. If there are any NaN
%       values in the dataset the algorithm will only take the first
%       non-NaN part.
%       E.g.: if data = [nan, 1, 2, 3, nan, nan, 4, nan, 5, 6, nan, 7] only
%       [1, 2, 3] will be used.
%   numberOfStates: Scalar number of different gaussians to fit.
%
% Return values:
%   mu: mean values of the gaussians
%   sigma: standard deviations of the gaussian.
%   numberOfIterations: number of iterations the algorithm performed.
%
%Calling schemes:
%   [mu, sigma] = GMM_CPP(data, numberOfStates);
%   [mu, sigma, numberOfIterations] = GMM_CPP(...);