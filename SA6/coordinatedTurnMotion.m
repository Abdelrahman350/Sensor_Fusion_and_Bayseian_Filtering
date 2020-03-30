function [fx, Fx] = coordinatedTurnMotion(x, T)
%COORDINATEDTURNMOTION calculates the predicted state using a coordinated
%turn motion model, and also calculated the motion model Jacobian
%
%Input:
%   x           [5 x 1] state vector
%   T           [1 x 1] Sampling time
%
%Output:
%   fx          [5 x 1] motion model evaluated at state x
%   Fx          [5 x 5] motion model Jacobian evaluated at state x
%
% NOTE: the motion model assumes that the state vector x consist of the
% following states:
%   px          X-position
%   py          Y-position
%   v           velocity
%   phi         heading
%   omega       turn-rate

% Your code here
px = x(1);
py = x(2);
v = x(3);
phi = x(4);
omega = x(5);

fx = [(px + T * v * cos(phi)), (py + T * v * sin(phi)), v, (phi + T * omega), omega]';

Fx = [1, 0, T*cos(phi), -T*v*sin(phi), 0;
      0, 1, T*sin(phi), T*v*cos(phi) , 0;
      0, 0, 1           , 0              , 0;
      0, 0, 0           , 1              , T;
      0, 0, 0           , 0              , 1];
end