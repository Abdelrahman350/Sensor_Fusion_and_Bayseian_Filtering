clc
clear
close all

tol = 1e-5;

for k=1:100
    % Random state dimension
    n = ceil(rand*10);
    % random mean
    x = rand(n,1);
    % random covariance
    P = rand(n,n);
    P = P*P';
    
    % UKF
    [SP,W] = sigmaPoints(x, P, 'UKF');
%    [~,W_ref] = reference.sigmaPoints(x, P, 'UKF');
        
    %Check the weights
%    assert(all(abs(W-W_ref)<tol),'Incorrect UKF weights')
    
    % CKF
    [SP,W] = sigmaPoints(x, P, 'CKF');
%    [~,W_ref] = reference.sigmaPoints(x, P, 'CKF');
        
    %Check the weights
    assert(all(abs(W-W_ref)<tol),'Incorrect CKF weights')
end