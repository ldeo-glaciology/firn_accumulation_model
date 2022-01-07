function dydz = odeSigW(zeta,y,n,m,Ar,r,beta)   
% this version of the ODE model has three simplifications compared to the original ODE model:
%  1) r(zeta) = r_s
% 2) sigma = -zeta
% 3) w = beta

% equation 25 in the manuscript 


% y(1) --> phi

dydz = -zeta^n*y(1)^m*(1-y(1))/Ar/beta/r;

end