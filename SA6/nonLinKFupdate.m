function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Predicted mean
%   P           [n x n] Predicted covariance
%   y           [m x 1] measurement vector
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state), 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%               Function must include all model parameters for the particular model, 
%               such as sensor position for some models.
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%
[hx,Hx]=h(x);
sqp = chol(P,'lower');

    switch type
        case 'EKF'
            
            % Your EKF update here
            S = Hx * P * Hx' + R;
            K = P * Hx' * S^(-1);
            P = P - K * S * K';
            x = x + K * (y - hx);

         case 'UKF'
    
            % Your UKF update here
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
            [SP,W] = sigmaPoints(x, P, type);
            yk = 0;
            Pxy = 0;
            Sk = 0;
            
            for i = 1:2*n+1
                yk = yk + W(i)*h(SP(:,i));
            end
            for i = 1:2*n+1
                Pxy = Pxy + (SP(:,i)-x) * (h(SP(:,i))-yk)' * W(i);
                Sk = Sk + (h(SP(:,i))-yk) * (h(SP(:,i))-yk)' * W(i);
            end
            Sk = R + Sk;
            x = x + Pxy * Sk^(-1) * (y - yk);
            P = P - Pxy * Sk^(-1) * Pxy';
        case 'CKF'
    
            % Your CKF update here
            [n,~] = size(x);
            k = 2*n;
            SP = zeros(n, k);
            W = zeros(1, 2*n);
            
            SP(:, 1) = x;
            W(1) = 1 - n/3;
            
            for i = 2:n
                W(i) = (1-W(1))/(2*n);
                W(i+n) = (1-W(1))/(2*n);
                SP(:,i) = x + sqrt(n / (1 - W(1)))*sqp(:,i-1);
                SP(:,i+n) = x - sqrt(n / (1 - W(1)))*sqp(:,i-1);
            end
            [SP,W] = sigmaPoints(x, P, type);
            yk = 0;
            Pxy = 0;
            Sk = 0;
            
            for i = 1:2*n
                yk = yk + W(i)*h(SP(:,i));
            end
            for i = 1:2*n
                Pxy = Pxy + (SP(:,i)-x) * (h(SP(:,i))-yk)' * W(i);
                Sk = Sk + (h(SP(:,i))-yk) * (h(SP(:,i))-yk)' * W(i);
            end
            Sk = R + Sk;
            x = x + Pxy * Sk^(-1) * (y - yk);
            P = P - Pxy * Sk^(-1) * Pxy';
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end

