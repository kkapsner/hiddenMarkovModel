%HMM_CPP performs a hidden markov model (HMM) analysis on the given dataset
%
%Parameter:
%   data: vector of the data to be modeled with a HMM
%   transitions: initial transition probability to change from state i to
%       state j: transition(i, j)
%       This has to be a MxM matrix. M defines the number of states that
%       are used in the model. The diagonal values are ignored and computed
%       internaly to have the sum of 1 in every column.
%   gaussianDefinitions: definition of the initial emission probabilities
%       by providing the mean and standard deviation for every state. This
%       has to be a Mx2 matrix where M matches the M of the transitions.
%       gaussianDefinitions(i, 1): mean of the gaussian for state i
%       gaussianDefinitions(i, 2): std of the gaussian for state i
%   options (optional): a struct that may have the following fields
%       (additional fields will be ignored and default values are in
%       brackets):
%       
%       verbose (true): if a verbose output should be displayed.
%       verboseOutputEmission (false): if the emission probabilities should
%           be included in the verbose output in every iteration.
%       verboseOutputTransition (true): if the transition probabilities
%           should be included in the verbose output in every iteration.
%       minSelfTransition (0): the minimal value the diagonal elements of
%           the transition probability matrix can have.
%       minEmision (1e-6): the minimal value an emission probability for
%           one bin can have.
%       doEmissionUpdate (true): if the emission probabilities should be
%           updated in every iteration step.
%       doTransitionUpdate (true): if the transition probabilities should
%           be updated in every iteration step.
%       binningCount (300): number of bins the data will be binned in.
%       maxiterations (100): maximum number of iterations
%       abortStateChanges (5): if less of equal data points change their
%           assigned state in the last iteration the algorithm will stop.
%       
%       options can also be a string containing a path to a JSON
%       configuration file. The file should look similar to:
%             {
%                 "verbose": {
%                     "enabled": false,
%                     "outputTransition": true,
%                     "outputEmission": false
%                 },
% 
%                 "minSelfTransition": 0,
%                 "minEmission": 1e-6,
% 
%                 "doEmissionUpdate": true,
%                 "doTransitionUpdate": true,
%                 "binningCount": 300,
%                 "maxIterations": 100,
%                 "abortStateChanges": 5
%             }
%
% Return values:
%   fittedStates: the fitted states to the data points
%   finalTransitions: final transition probabilities after the last
%       iteration. This will be a MxM matrix where M is the number of
%       states in the model.
%   finalEmissions: final emission probabilities after the last iteration.
%       This will be a MxN matrix where M is the number of states in the
%       model and N is the number of bins where the data was binned in (see
%       binningCount field in options parameter).
%   numberOfIterations: number of iterations the algorithm performed.
%
%Calling schemes:
%   fittedStates = HMM_CPP(data, transitions, gaussianDefinitions);
%   fittedStates = HMM_CPP(..., options);
%   [fittedStates, finalTransitions, finalEmissions, numberOfIterations] =
%       HMM_CPP(...);
%
% example code:
%
% means = [1; 4; 7]; % make sure this is a column vector
% std = [0.1; 0.2; 1]; % this also has to be a column vector
% transitions = [ ...
%     [0; 0.1; 0.1], ... the 0 will be computed to 0.8
%     [0.1; 0; 0], ... the first 0 will be computed to 0.9
%     [0.3; 0; 0] ... the seconde 0 will be computed to 0.7
% ];
% 
% options = struct( ...
%     'verbose', true, ...
%     'verboseOutputEmission', false, ...
%     'verboseOutputTransition', true, ...
%     'minSelfTransition', 0, ...
%     'minEmision', 1e-6, ...
%     'doEmissionUpdate', true, ...
%     'doTransitionUpdate', true, ...
%     'binningCount', 300, ...
%     'maxiterations', 100, ...
%     'abortStateChanges', 5 ...
% );
% [states, transitionMatrix, emissionProbabilties, iterationCount] = ...
%     HMM_cpp(data, transtions, [means, stds], options);

