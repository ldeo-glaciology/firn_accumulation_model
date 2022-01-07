function dydz = odeSig_r(zeta,y,n,m,Ar)   
% this version of the ODE model has two simplifications compared to the original ODE model:
% 1) delta = 0
% 2) sigma = -z

% equation 26 in the manuscript 

% y(1) --> phi
% y(2) --> w
% y(3) --> r

dydz = [-zeta^n*y(1)^m*(1-y(1))/Ar/y(2)/y(3)
    -zeta^n*y(1)^m/Ar/y(3)
    1/y(2)];
end