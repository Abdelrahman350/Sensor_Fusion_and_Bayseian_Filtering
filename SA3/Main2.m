clc
clear
close all

absTol = 1e-1;
relTol = 5e-2;


N = 10000;

n = randi(5,1);
m = randi(n,1);

N = m*N;

% Define state sequence
X = rand(n,N+1);

% Define measurement model
H = randn(m,n);
V = rand(m,m);
R = V*diag(10*rand(m,1))*V';

% Generate measurements
Y = genLinearMeasurementSequence(X, H, R);

assert(size(Y,1) == m, 'Y has the wrong measurement dimension');
assert(size(Y,2) == N, 'Y should have N columns');

Rest = cov((Y-H*X(:, 2:N+1))');
assert(all(all((mean(Y-H*X(:, 2:N+1),2) < absTol))), 'Measurement noise is not zeros mean');
assert(all(all((abs(Rest-R) < relTol*R))), 'Measurement noise covariance is not within tolerances');
