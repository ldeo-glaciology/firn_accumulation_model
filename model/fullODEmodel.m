function dydz = fullODEmodel(zeta,y,Ar,delta)   % full ODE model

%% encodes the ODE steady state firn model (section 2.4 of the manuscript)

% example of how to solve it:
%     options = odeset('RelTol',1e-10,'AbsTol',1e-10);
%     [zeta,ODEModel] = ode45(@(x,y) fullODEmodel(x,y,Ar,delta),[0 1],[phi_s 0 -beta/(1-phi_s) r_s 0],options);
%

%% variable order 
% y(1) phi
% y(2) sigma
% y(3) w
% y(4) r
% y(5) A
dydz = [abs(y(2))*y(1)*(1-y(1))/(Ar*y(3)*y(4))
    -(1-y(1))
    (abs(y(2))*y(1))/(Ar*y(4))
    -(1-delta*y(4))/y(3)
    -1/y(3)];
end
