function [mu, sigma2] = posteriorGaussian_n(mu_x, sigma2_x, y, sigma2_r)
%posteriorGaussian performs a single scalar measurement update with a
%measurement model which is simply "y = x + noise".
%
%Input
%   MU_P            The mean of the (Gaussian) prior density.
%   SIGMA2_P        The variance of the (Gaussian) prior density.
%   SIGMA2_R        The variance of the measurement noise.
%   Y               The given measurement.
%
%Output
%   MU              The mean of the (Gaussian) posterior distribution
%   SIGMA2          The variance of the (Gaussian) posterior distribution
n = 16;
%Your code here
mu = sigma2_r * ( mu_x / sigma2_x + n * y / sigma2_r);
sigma2 = (1./sigma2_x + n ./sigma2_r).^-1;

end