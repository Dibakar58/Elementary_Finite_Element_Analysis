function [shape,nderiv] = shape2d(xi,eta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
shape=(1/4)*[(1-xi)*(1-eta);
    (1+xi)*(1-eta);
    (1+xi)*(1+eta);
    (1-xi)*(1+eta)] ;

nderiv=(1/4)*[-(1-eta),-(1-xi);
    (1-eta),-(1+xi);
    (1+eta),(1+xi);
    -(1+eta),(1-xi) ] ;


end

