function X = genLinearStateSequence(x_0, P_0, A, Q, N)
%GENLINEARSTATESEQUENCE generates an N-long sequence of states using a 
%    Gaussian prior and a linear Gaussian process model
%
%Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   N           [1 x 1] Number of states to generate
%
%Output:
%   X           [n x N+1] State vector sequence
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Your code here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x0 = mvnrnd(xk, P_0)';
%qk = zeros(size(x_0));
[n, ~]=size(x_0);
X = zeros(n, N+1);
X(:,1) = mvnrnd(x_0, P_0)';
for k = 2: N+1
    X(:,k) = A*X(:,k-1) + mvnrnd(zeros(size(x_0)), Q)';
end
end