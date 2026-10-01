function [MP1,CMC1] = mpath(C1C, L1C, f1, L2W, f2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

c = 299792458; %m/s

lambda1 = c./f1;
lambda2 = c./f2;

MP1 = C1C - ((f1^2 + f2^2)./(f1^2 - f2^2)).*L1C.*lambda1  + ((2*f2^2)./(f1^2 - f2^2)).*L2W.*lambda2;

CMC1 = C1C - L1C.*lambda1;

end 

