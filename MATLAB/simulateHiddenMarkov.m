function [emission, states] = simulateHiddenMarkov(transition, mean, std, dataSize)
    statesCount = numel(mean);
    emission = zeros(dataSize, 1);
    states = zeros(dataSize, 1);
    transition = cumsum(transition, 2);
    
    states(1) = randi(statesCount);
    for i = 2:dataSize
        r = rand(1);
        states(i) = find(r < transition(states(i-1), :), 1, 'first');
    end
    for i = 1:dataSize
       emission(i) = normrnd(mean(states(i)), std(states(i)));
    end
end

