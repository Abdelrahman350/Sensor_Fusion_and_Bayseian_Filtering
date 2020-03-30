function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, ...
                             N, bResample, plotFunc)
%PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
% state-space model.
%
% Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   Y           [m x K] Measurement sequence to be filtered
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%   N           Number of particles
%   bResample   boolean false - no resampling, true - resampling
%   plotFunc    Handle for plot function that is called when a filter
%               recursion has finished.
% Output:
%   xfp         [n x K] Posterior means of particle filter
%   Pfp         [n x n x K] Posterior error covariances of particle filter
%   Xp          [n x N x K] Particles for posterior state distribution in times 1:K
%   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K

% Your code here, please. 
% If you want to be a bit fancy, then only store and output the particles if the function
% is called with more than 2 output arguments.
[m, K] = size(Y);
[n, ~] = size(x_0);
xfp = zeros(n, K);
Pfp = zeros(n, n, K);
Xp  = zeros(n, N, K);
Wp  = zeros(N, K);

X_r = mvnrnd(x_0, P_0)' .* ones(n,N);
Wr = 1/N .* ones(1,N);

for i = 1:K
    [X_k, W_k] = pfFilterStep(X_r, Wr, Y(:,i), proc_f, proc_Q, meas_h, meas_R);
    [Xr, Wr, j] = resampl(X_k, W_k);
    mean_nRS = Xr * Wr';%sum(Xr .* Wr,2);
    cov_nRS = (Xr - mean_nRS) * ((Xr - mean_nRS)' .* Wr');
    xfp(:,i) = mean_nRS;
    Pfp(:,:,i) = cov_nRS;
    Xp(:,:,i) = Xr;
    Wp(:,i) = Wr';
end
plotFunc;
end

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
[~,ind1] = sort([u,wc]);
ind2 = find(ind1<=N);
j = ind2 - (0:N-1);
Xr = X(:,j);
Wr = ones(1,N)./N;
end

function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)
%PFFILTERSTEP Compute one filter step of a SIS/SIR particle filter.
%
% Input:
%   X_kmin1     [n x N] Particles for state x in time k-1
%   W_kmin1     [1 x N] Weights for state x in time k-1
%   y_k         [m x 1] Measurement vector for time k
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%
% Output:
%   X_k         [n x N] Particles for state x in time k
%   W_k         [1 x N] Weights for state x in time k

% Your code here!
[n, N] = size(X_kmin1);
[~, m] = size(meas_R);
X_k = zeros(n,N);
W_k = zeros(1,N);
py_x = zeros(1,N);
%prediction

X_k =  proc_f(X_kmin1) + mvnrnd(zeros(n,1), proc_Q, N)';

for k = 1: N
    %X_k(:,k) =  proc_f(X_kmin1(:,k)) + mvnrnd(zeros(n,1), proc_Q)';
    py_x(:,k) = exp(-0.5*(yk-meas_h(X_k(:,k)))' * meas_R^(-1) * (yk-meas_h(X_k(:,k))));
end
W_k = W_kmin1 .* py_x;
W_k = W_k./ sum(W_k);
end