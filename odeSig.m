function dydz = odeSig(zeta,y,n,m,Ar,r)   
% this version of the ODE model has two simplifications, compared to the original ODE model:
% 1) r(zeta) = r_s
% 2) sigma = zeta

% example of how to run this
% [zeta,y_ODE1] = ode45(@(x,y) odeSig(x,y,n,m,Ar,r_s), p.z_h,[phi_s beta(jj)/(1-phi_s)],options);


% equation 24 in the manuscript 

% y(1) --> phi
% y(2) --> w

dydz = [-zeta^n*y(1)^m*(1-y(1))/Ar/y(2)/r
    -zeta^n*y(1)^m/Ar/r];
end