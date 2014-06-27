function states = HMM(data, dim, transitions, gaussianDefinitions, options)
%HMM is a wrapper for HMM_CPP that allows a multidimensional input
%
%   states = HMM(data, transitions, gaussianDefinitions)
%   states = HMM(data, DIM, transitions, gaussianDefinitions)
%   states = HMM(..., options)
%
%The hidden Markov analysis is only performed along one dimension of the
%data matrix. If no DIM parameter is provided the first non singular
%dimension is used. Otherwise the HMM runs along the dimension specified in
%DIM. All other parameter are the same as for HMM_CPP
%
%SEE ALSO: HMM_CPP
    
    if (nargin < 5)
        if (size(transitions, 1) ~= size(gaussianDefinitions, 1))
            if (nargin == 4)
                options = gaussianDefinitions;
            else
                options = [];
            end
            if (nargin >= 3)
                gaussianDefinitions = transitions;
                transitions = dim;
                dim = [];
            end
        else
            options = [];
        end
    end
    
    if (isempty(dim))
        [data, nshifts] = shiftdim(data);
    elseif (dim ~= 1)
        ndim = ndims(data);
        assert(isscalar(dim) && dim > 0 && dim <= ndim, ...
            'HMM:InvalidDimensions');
        perm = [dim, 1:(dim-1), (dim+1):ndim];
        data = permute(data, perm);
    end
    
%     states = zeros(size(data));
%     numBlocks = numel(data) / size(data, 1);
%     
%     for i = 1:numBlocks
%         states(:, i) = HMM_cpp(data(:, i), transitions, gaussianDefinitions, options);
%     end
    states = reshape( ...
        HMM_cpp(data, transitions, gaussianDefinitions, options), ...
        size(data) ...
    );
    
    
    if (isempty(dim))
        states = shiftdim( ...
            states, ...
            -nshifts ...
        );
    elseif (dim ~= 1)
        states = ipermute( ...
            states, ...
            perm ...
        );
    end 
end