function dydz = odeSigW_r(zeta,y,n,m,Ar,r_s,beta)   
% this version of the ODE model was three simplifications compared to the original ODE model:
% 1) delta = 0
% 2) sigma = -z
% 3) w = beta

% equation 27 in the manuscript 


% y(1) --> phi

dydz = -zeta^(n-1)*y(1)^m*(1-y(1))/Ar/(1+(beta*r_s/zeta));
end