clc
clear
close all

tol = 1e-5;


n = 1;
T = 1;

% Motion Parameters
A = 2;
Q = T*1;

% Prior
xPrior  = mvnrnd(zeros(n,1)', diag(ones(n)))';
V       = rand(n,n);
PPrior  = V*diag(rand(n,1))*V';

[xPred, PPred] = linearPrediction(xPrior, PPrior, A, Q);

figure(1); clf; hold on;
x = linspace(xPred - 4*sqrt(PPred), xPred + 4*sqrt(PPred),100);
plot(x, normpdf(x,xPred, sqrt(PPred)));
plot(x, normpdf(x,xPrior, sqrt(PPrior)));
title('Your solution')
xlabel('x');
ylabel('p(x)')
legend('Predicted density', 'Prior density');

figure(2); clf; hold on;
x = linspace(xPred_ref - 4*sqrt(PPred_ref), xPred_ref + 4*sqrt(PPred_ref),100);
plot(x, normpdf(x, xPred_ref, sqrt(PPred_ref)));
plot(x, normpdf(x, xPrior, sqrt(PPrior)));
title('Reference solution')
xlabel('x');
ylabel('p(x)')
legend('Predicted density', 'Prior density');


assert(isequal(size(xPrior),size(xPred)), 'Dimenstion of prior and predicted mean need to be the same.');
assert(isequal(size(PPrior),size(PPred)), 'Dimenstion of prior and predicted covariance need to be the same.');
assert(all(abs(xPred-xPred_ref)<tol), 'Predicted mean is not within tolerance.');
assert(all(all(abs(PPred-PPred_ref)<tol)), 'Predicted covarinace is not within tolerance.');

[~, p] = chol(PPred);
assert(p == 0 || trace((PPred)) ~= 0, 'Posterior covariance is not positive semi definite covarinace matrix');
