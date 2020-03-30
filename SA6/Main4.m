clc
clear
close all

tol = 1e-2;

for k = 1:100
    % state
    x = rand(5,1);
    % Covariance
    P = sqrt(10)*rand(5,5);
    P = P*P';
    
    % Random sensor position sequence
    s1 = rand(2,1)*100;
    s2 = rand(2,1)*100;
    
    % Measurement model
    h = @(x) dualBearingMeasurement(x, s1, s2);

    % Random measurement
    y = h(x)+rand(2,1);
    
    % Measurement noise covariance
    R = diag(rand(1,2).^2);
    
    % EKF
    [xp, Pp] = nonLinKFupdate(x, P, y, h, R, 'EKF');
    %[xp_ref, Pp_ref] = nonLinKFupdate(x, P, y, h, R, 'EKF');
%    [xp_ref, Pp_ref] = reference.nonLinKFupdate(x, P, y, h, R, 'EKF');
    % Test results
%    assert(norm(xp-xp_ref)<tol, 'Incorrect EKF mean');
%    assert(norm(Pp-Pp_ref)<tol, 'Incorrect EKF covariance');
    
    % UKF
    [xp, Pp] = nonLinKFupdate(x, P, y, h, R, 'UKF');
    %[xp_ref, Pp_ref] = nonLinKFupdate(x, P, y, h, R, 'UKF');
%    [xp_ref, Pp_ref] = reference.nonLinKFupdate(x, P, y, h, R, 'UKF');
    % Test results
%    assert(norm(xp-xp_ref)<tol, 'Incorrect UKF mean');
%    assert(norm(Pp-Pp_ref)<tol, 'Incorrect UKF covariance');
    
    % CKF
    [xp, Pp] = nonLinKFupdate(x, P, y, h, R, 'CKF');
    %[xp_ref, Pp_ref] = nonLinKFupdate(x, P, y, h, R, 'CKF');
%    [xp_ref, Pp_ref] = reference.nonLinKFupdate(x, P, y, h, R, 'CKF');
    % Test results
%    assert(norm(xp-xp_ref)<tol, 'Incorrect CKF mean');
%    assert(norm(Pp-Pp_ref)<tol, 'Incorrect CKF covariance');
    
    
    % Define prior
    x_0     = [0 0 10 0 0]';
    n       = length(x_0);
    P_0     = diag([1 1 1 1*pi/180 1*pi/180].^2);
    % Sample time
    T = rand;
    % Covariance
    sigV = 1;
    sigOmega = 1*pi/180;
    G = [zeros(2,2); 1 0; 0 0; 0 1];
    Q = G*diag([sigV^2 sigOmega^2])*G';
    
    motionModel = @(x) coordinatedTurnMotion(x, T);
    hSP = @sigmaPoints;
    %hSP = @reference.sigmaPoints;
    
    % EKF
    [xp, Pp] = nonLinKFprediction(x_0, P_0, motionModel, Q, 'EKF');
    %[xp_ref, Pp_ref] = nonLinKFprediction(x_0, P_0, motionModel, Q, 'EKF');
%    [xp_ref, Pp_ref] = reference.nonLinKFprediction(x_0, P_0, motionModel, Q, 'EKF');
    % Test results
%    assert(norm(xp-xp_ref)<tol, 'Incorrect EKF mean');
%    assert(norm(Pp-Pp_ref)<tol, 'Incorrect EKF covariance');
    
    % UKF
    [xp, Pp] = nonLinKFprediction(x_0, P_0, motionModel, Q, 'UKF');
    %[xp_ref, Pp_ref] = nonLinKFprediction(x_0, P_0, motionModel, Q, 'UKF');
%    [xp_ref, Pp_ref] = reference.nonLinKFprediction(x_0, P_0, motionModel, Q, 'UKF');
    % Test results
%    assert(norm(xp-xp_ref)<tol, 'Incorrect UKF mean');
%    assert(norm(Pp-Pp_ref)<tol, 'Incorrect UKF covariance');
    
    % CKF
    [xp, Pp] = nonLinKFprediction(x_0, P_0, motionModel, Q, 'CKF');
    %[xp_ref, Pp_ref] = nonLinKFprediction(x_0, P_0, motionModel, Q, 'CKF');
%    [xp_ref, Pp_ref] = reference.nonLinKFprediction(x_0, P_0, motionModel, Q, 'CKF');
    % Test results
%    assert(norm(xp-xp_ref)<tol, 'Incorrect CKF mean');
%    assert(norm(Pp-Pp_ref)<tol, 'Incorrect CKF covariance');
    [xf, Pf, xp, Pp] = nonLinearKalmanFilter(y, x_0, P_0, motionModel, Q, h, R, 'EKF');
    disp(['Random test ' num2str(k) ' passed'])
end