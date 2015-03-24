function [l,lo]=lineParam(A,gammaR)

%lineParam - Create the equation for a line in space
%
%
%[l,lo]=lineParam(A,alphaR)
%       Inputs:
%               A - Tranformation matrix (4x4)
%               gammaR - length of line
%
%       Outputs:
%               l - equation of a line in space
%

%Create a basis starting vector
p=[1 0 0]';

%Scale it by length
lo=gammaR*p;

%Rotate the line
l=A*[lo;1];

%Upper 3 contain the end point equation
l=l(1:3);