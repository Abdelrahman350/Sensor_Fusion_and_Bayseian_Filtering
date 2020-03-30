function [mu_y, Sigma_y, y_s] = approxGaussianTransform(mu_x, Sigma_x,f, N)
%approxGaussianTransform takes a Gaussian density and a transformation 
%function and calculates the mean and covariance of the transformed density.
%
%Inputs
%   MU_X        [m x 1] Expected value of x.
%   SIGMA_X     [m x m] Covariance of x.
%   F           [Function handle] Function which maps a [m x 1] dimensional
%               vector into another vector of size [n x 1].
%   N           Number of samples to draw. Default = 5000.
%
%Output
%   MU_Y        [n x 1] Approximated mean of y.
%   SIGMA_Y     [n x n] Approximated covariance of y.
%   ys          [n x N] Samples propagated through f


if nargin < 4
    N = 5000;
end

%Your code here
[m,n] = size(mu_x);
%x = zeros(m,N);
%for i = 1:m
    %x(i,:) = linspace(-3*Sigma_x(i)+mu_x(i),3*Sigma_x(i)+mu_x(i),N);
%end
%X = zeros(m,N);

%for i = 1:N
    %X(:,i) = 1/sqrt(det(2*pi*Sigma_x))*exp(-0.5*(x(:,i)-mu_x)'*Sigma_x^(-1)*(x(:,i)-mu_x));
%end
X = mvnrnd(mu_x,Sigma_x,N)';

figure
plot(X(1,:),X(2,:),'+')
f = @func2;
y_s = f(X);
mu_y = sum(y_s,2)/N;
[my,ny]=size(mu_y);
S = zeros(my,my);
for i = 1:N
    S = S + (y_s(:,i)-mu_y)*(y_s(:,i)-mu_y)';
end

Sigma_y = S/(N-1);
end
function f = func2(X)
    x=X(1,:);
    y=X(2,:);
    f(1,:) = sqrt(x.^2 + y.^2);
    f(2,:) = atan2(y,x);
end
