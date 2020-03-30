function [X, P] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x N] Estimated state vector sequence
%   P           [n x n x N] Filter error convariance
%

%% Parameters
N = size(Y,2);

n = length(x_0);
m = size(Y,1);

%% Data allocation
X = zeros(n,N);
P = zeros(n,n,N);
X(:,1) = mvnrnd(x_0, P_0)';
P(:,:,1) = P_0;
[X(:,1), P(:,:,1)] = linearPrediction(x_0, P_0, A, Q);
[X(:,1), P(:,:,1)] = linearUpdate(X(:,1), P(:,:,1), Y(:,1), H, R);
for k = 2: N
    
    [X(:,k), P(:,:,k)] = linearPrediction(X(:,k-1), P(:,:,k-1), A, Q);
    
    [X(:,k), P(:,:,k)] = linearUpdate(X(:,k), P(:,:,k), Y(:,k), H, R);
    %X(:,k) = A*X(:,k-1) + mvnrnd(zeros(size(x_0)), Q)';
end
   
end

