function [x, P] = nonLinKFprediction(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%
sqp = chol(P,'lower');
    switch type
        case 'EKF'
            
            % Your EKF code here
            [fx, Fx] = f(x);
            x = fx;
            P = Q + Fx * P * Fx'; 
            
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
            
            x = 0;
            P = 0;
            for i = 1:2*n+1
                x = x + W(i)*f(SP(:,i));
            end
            for i = 1:2*n+1
                P = P + (f(SP(:,i))-x) * (f(SP(:,i))-x)' * W(i);
            end
            P = P + Q;
                        
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
            
            x = 0;
            P = 0;
            for i = 1:2*n
                x = x + W(i)*f(SP(:,i));
            end
            for i = 1:2*n
                P = P + (f(SP(:,i))-x) * (f(SP(:,i))-x)' * W(i);
            end
            P = P + Q;
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end