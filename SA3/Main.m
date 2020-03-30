clc
clear
close all

% Tolerance
tol = 1e-1;

N = 5000;

% Define prior
x_0     = [0;0]; 
n       = length(x_0); 
P_0     = diag(ones(n,1));

% Define process model
A       = [1 1; 0 1];
Q       = diag(ones(n,1));

% generate state sequence
s = rng;
X = genLinearStateSequence(x_0, P_0, A, Q, N);
rng(s);

% Whiteness
Qest = cov((X(:,2:end)-A*X(:,1:end-1))');

% Plot results
figure(2);clf;hold on;
subplot(2,1,1);plot(X(1,:));
title('Your solution');
xlabel('k');
ylabel('x-position');
subplot(2,1,2);plot(X(2,:));
xlabel('k');
ylabel('speed');