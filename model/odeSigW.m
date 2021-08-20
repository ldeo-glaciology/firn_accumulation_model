function dydz = odeSigW(zeta,y,n,m,Ar,r,beta)   
% this version of the ODE model has three simplifications compared to the original ODE model:
% 1) delta = 0
% 2) sigma = -zeta
% 3) w = -beta

% equation 26 in the manuscript (check this later)



% y(1) --> phi

dydz = -zeta^n*y(1)^m*(1-y(1))/Ar/beta/r;

end