function [l,lo]=lineParam(A5_0,alphaR)

p=[1 0 0]';
lo=alphaR*p;
l=A5_0*[lo;1];
l=l(1:3);