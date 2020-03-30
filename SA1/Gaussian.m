clc
clear
close all

mu_x=[7;3];
Sigma_x=[0.2 0;0 8];
f = @func2;
[mu_y, Sigma_y, ~] = approxGaussianTransform(mu_x, Sigma_x,f,2000)