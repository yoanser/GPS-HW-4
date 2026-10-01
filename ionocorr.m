function [PRIF,iono] = ionocorr(firstPseudo,f1,secondPseudo,f2)
%{
Author: Yoan Serafimov
INPUTS:
firstPseudo - first pseudorange measurment columnn vector m
f1 - frequency associated with first pseudo Hz
secondPseudo - second pseudorange column vector m
f2 - frequency associated with second pseudo Hz

OUTPUTS:
PRIF = ionofree pseudorange (m)
iono = iono error (m)
%}
PRIF = (f1^2./(f1^2-f2^2)).*firstPseudo  -  (f2^2./(f1^2-f2^2)).*secondPseudo;
iono = (f2^2./(f1^2-f2^2)).*(secondPseudo - firstPseudo);
end

