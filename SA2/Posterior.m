clc
clear

y = [200, 195];
sigma2_r = 40^2;
sigma2_x = [10^2, 60^2];
mu_x = [180, 160];
[mu, sigma2] = posteriorGaussian_n(mu_x, sigma2_x, y, sigma2_r)
sqrt(sigma2)