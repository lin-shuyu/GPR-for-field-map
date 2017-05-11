%% Create a synthetic model of the variation in RSSi of one AP (router) over
% a space of 20*20 (m)

%% inputs:
% x - any input position as specified, 2*1 column vector (horiztal;vertical)
% x_AP - the position of the AP in the room, 2*1 column vector (horiztal;vertical)
% n - Gaussian white noise @ x

%% outputs
% f - the noise-free function value (RSSi in this case, in dBm) @ x
% y - white gaussian noise added to f

function [f,y] = RSSi(x,x_AP,e)
% define a few constants in the formula
A = -30; % the received signal strength @ 1m
n = 3; % path-loss exponent

d = norm(x - x_AP);
f = -10*n*log10(d) + A; 
y = f + e*randn;
end