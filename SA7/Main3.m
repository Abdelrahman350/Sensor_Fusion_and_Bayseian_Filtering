clc
clear
close all
% Test PF filter and compare with Kalman filter output. Generate 2D state sequence (CV, pos and vel) and 1D measurement sequence (pos) and compare
% outputs with Kalman filter.

% Set prior
sigma = 2;
x_0 = [0 1]';
P_0 = [sigma^2 0; ...
       0 sigma^2];
n = size(x_0,1);

% Number of time steps
K = 20;

% Models
A = [1 0.1; 0 1];
Q = [0 0; 0 0.5];
H = [1 0];
R = 1;

m = 1;

% Generate state and measurement sequences
X = zeros(n,K);
Y = zeros(m,K);
   
q = mvnrnd([0 0], Q, K)';
r = mvnrnd(zeros(1,m), R, K)';
x_kmin1 = x_0;
for k = 1:K
    xk = f(x_kmin1,A) + q(:,k);
    X(:,k) = xk;
    x_kmin1 = xk;
    
    Y(:,k) = h(xk, H) + r(:,k);
end

% Run Kalman filter
[Xf, Pf] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

% Run PF filter with and without resampling
N = 20000;
proc_f = @(X_kmin1) (f(X_kmin1, A));
meas_h = @(X_k) (h(X_k, H));
plotFunc = @(k, Xk, Xkmin1, Wk, j) (0); % Define dummy function that does nothing

[xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, Q, meas_h, R, ...
                              N, false, plotFunc);
                          
[xfpr, Pfpr, Xpr, Wpr] = pfFilter(x_0, P_0, Y, proc_f, Q, meas_h, R, ...
                                  N, true, plotFunc);


plot(X(2,:));
plot(X(1,:));
plot(X(1,:), X(2,:));
hold on;
plot(Xf(2,:), 'c');
plot(Xf(1,:), 'c');
plot(Xf(1,:), Xf(2,:), 'c');
plot(xfp(2,:), 'r');
plot(xfpr(2,:), 'g');
plot(xfp(1,:), 'r');
plot(xfpr(1,:), 'g');
plot(xfp(1,:), xfp(2,:), 'r');
plot(xfpr(1,:), xfpr(2,:),'g');
hold off;

% Compute means and covariances for each particle filter output

for k = 1:K
    mean_nRS = Xp(:,:,k) * Wp(:,k);
    cov_nRS = (Xp(:,:,k) - mean_nRS) * ((Xp(:,:,k) - mean_nRS)' .* Wp(:,k));
    % Compare with the provided means and covariances
    assert(norm(mean_nRS-xfp(:,k)) < 0.00001, ...
        'The output mean from PF without RS is not consistent with the output particles.');
    assert(norm(cov_nRS-Pfp(:,:,k)) < 0.00001, ...
        'The output covariance from PF without RS is not consistent with the output particles.');
    
    mean_RS = Xpr(:,:,k) * Wpr(:,k);
    cov_RS = (Xpr(:,:,k) - mean_RS) * ((Xpr(:,:,k) - mean_RS)' .* Wpr(:,k));
    % Compare with the provided means and covariances
    assert(norm(mean_RS-xfpr(:,k)) < 0.00001, ...
        'The output mean from PF with RS is not consistent with the output particles.');
    assert(norm(cov_RS-Pfpr(:,:,k)) < 0.00001, ...
        'The output covariance from PF with RS is not consistent with the output particles.');

    % Compare with Kalman filter
%     assert(norm(Xf(:,k)-xfp(:,k)) < 0.5, ...
%         'The mean from PF WITHOUT resampling deviates too much from KF filter mean.');  
%     assert(norm(Pf(:,:,k)-Pfp(:,:,k)) < 1, ...
%         'The covariance from PF WITHOUT resampling deviates too much from KF filter covariance.');  

    assert(norm(Xf(:,k)-xfpr(:,k)) < 0.5, ...
        'The mean from PF WITH resampling deviates too much from KF filter mean.');  
    assert(norm(Pf(:,:,k)-Pfpr(:,:,k)) < 1, ...
        'The covariance from PF WITH resampling deviates too much from KF filter covariance.'); 

end




function X_k = f(X_kmin1, A)
%
% X_kmin1:  [n x N] N states vectors at time k-1
% A:        [n x n] Matrix such that x_k = A*x_k-1 + q_k-1
    X_k = A*X_kmin1;
end

function H_k = h(X_k, H)
%
% X_k:  [n x N] N states
% H:    [m x n] Matrix such that y = H*x + r_k
    H_k = H*X_k;
end

