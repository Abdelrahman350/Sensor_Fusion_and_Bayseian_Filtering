function [hx, Hx] = dualBearingMeasurement(x, s1, s2)
%DUOBEARINGMEASUREMENT calculates the bearings from two sensors, located in 
%s1 and s2, to the position given by the state vector x. Also returns the
%Jacobian of the model at x.
%
%Input:
%   x           [n x 1] State vector, the two first element are 2D position
%   s1          [2 x 1] Sensor position (2D) for sensor 1
%   s2          [2 x 1] Sensor position (2D) for sensor 2
%
%Output:
%   hx          [2 x 1] measurement vector
%   Hx          [2 x n] measurement model Jacobian
%
% NOTE: the measurement model assumes that in the state vector x, the first
% two states are X-position and Y-position.

% Your code here
[n,~] = size(x);
Hx = zeros(2,n);

xk = x(1);
yk = x(2);

x1 = s1(1);
y1 = s1(2);

x2 = s2(1);
y2 = s2(2);

hx = [atan2(yk-y1,xk-x1); 
      atan2(yk-y2,xk-x2)];

Hx = [-(yk-y1)/((yk-y1)^2+(xk-x1)^2), (xk-x1)/((yk-y1)^2+(xk-x1)^2), 0, 0;
    -(yk-y2)/((yk-y2)^2+(xk-x2)^2), (xk-x2)/((yk-y2)^2+(xk-x2)^2), 0, 0];
end