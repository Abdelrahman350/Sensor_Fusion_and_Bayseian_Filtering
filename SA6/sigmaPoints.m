function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%
sqp = chol(P,'lower');
    switch type        
        case 'UKF'
    
            % your code
            [n,~] = size(x);
            k = 2*n + 1;
            SP = zeros(n, k);
            W = zeros(1, 2*n+1);
            
            SP(:, 1) = x;
            W(1) = 1 - n/3;
            
            for i = 2:n+1
                W(i) = (1-W(1))/(2*n);
                W(i+n) = (1-W(1))/(2*n);
                SP(:,i) = x + sqrt(n / (1 - W(1)))*sqp(:,i-1);
                SP(:,i+n) = x - sqrt(n / (1 - W(1)))*sqp(:,i-1);
            end
        case 'CKF'
            
            % your code
            [n,~] = size(x);
            k = 2*n;
            SP = zeros(n, k);
            W = zeros(1, 2*n);
            
            for i = 1:n
                W(i) = (1)/(2*n);
                W(i+n) = (1)/(2*n);
                SP(:,i) = x + sqrt(n)*sqp(:,i);
                SP(:,i+n) = x - sqrt(n)*sqp(:,i);
            end
        otherwise
            error('Incorrect type of sigma point')
    end

end