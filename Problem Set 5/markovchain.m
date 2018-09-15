function [chain] = markovchain(prob, T, start)
%markovchain Generate Markov chain
%   [chain] = markovchain(prob, T, start) returns chain, Markov chain.

chain = ones(T,1);  % Initialize Markov Chain
chain(1)= start;    % Set starting value

cum_prob = cumsum(prob,2);  % Compute cumulative distribution

% Generate Markov Chain using random numbers uniformly distributed
for t = 2:T
    chain(t)=find(cum_prob(chain(t-1),:)>rand(),1);
end

end