function [Xr, Wr, j] = resampl(X, W)
%RESAMPLE Resample particles and output new particles and weights.
% resampled particles. 
%
%   if old particle vector is x, new particles x_new is computed as x(:,j)
%
% Input:
%   X   [n x N] Particles, each column is a particle.
%   W   [1 x N] Weights, corresponding to the samples
%
% Output:
%   Xr  [n x N] Resampled particles, each corresponding to some particle 
%               from old weights.
%   Wr  [1 x N] New weights for the resampled particles.
%   j   [1 x N] vector of indices refering to vector of old particles

% Your code here!
N = length(W);
u = ([0:N-1]+rand(1))/N;
wc = cumsum(W);
wc = wc/wc(N);
[dum,ind1] = sort([u,wc]);
ind2 = find(ind1<=N);
j = ind2 - (0:N-1);
Xr = X(:,j);
Wr = ones(1,N)./N;
end
