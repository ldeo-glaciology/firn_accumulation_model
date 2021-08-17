function dydz = odeSig(zeta,y,n,m,Ar,r)   
% this version of the ODE model has two simplifications, compared to the original ODE model:
% 1) r(zeta) = r_s
% 2) sigma = -zeta

% equation 25 in the manuscript (check this later)

% y(1) --> phi
% y(2) --> w

dydz = [zeta^n*y(1)^m*(1-y(1))/Ar/y(2)/r
    zeta^n*y(1)^m/Ar/r];
end