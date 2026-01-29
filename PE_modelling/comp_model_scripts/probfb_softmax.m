function [probs] = probfb_softmax(Q,beta)

% INPUT:
% Q is the two estimated action values, e.g. [0.5 0.5]
% beta is the inverse temperature
% OUTPUT:
% the probability of the choices

probs = (exp(beta*Q))/(exp(beta*Q(1))+exp(beta*Q(2)));