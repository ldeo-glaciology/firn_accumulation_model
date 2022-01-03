function Out = FCM9b(p)

N = p.GridNumber;
D1 = p.Gradient;
z_h = p.z_h;
beta = p.beta;
xi_p = p.xi_p;
dxi_p = p.dxi_p;
delta = p.delta;
lambda_c = p.lambda_c;
lambda_g = p.lambda_g;
m = p.PorosityExponent;
n = p.StressExponent;
Pe = p.PecletNumber;
Fl = p.FluxNumber;
Ar = p.ArthenNumber;
w_s = p.w_s;
h_0 = p.h_0;
Vars_init = [p.InitialConditions(1:2*N); p.InitialConditions(3*N+1:end)]; 
sim_r = p.sim_r;
beta_multiplier = p.beta_multiplier;

%%% Solve the system
time = [0 10];
options = odeset('Events', @SteadyState, 'RelTol',1e-8,'AbsTol', 1e-10);
[Time, Vars] = ode15s(@GovEq, time, Vars_init, options);

%%% Collect the simulation variables.
Phi = Vars(:, 1:N)';
Rho = p.rho_i*(1 - Phi);
GrainSize = Vars(:, N+1:2*N)';
Age = Vars(:, 2*N+1:3*N)';
Temp = Vars(:, 3*N+1:4*N)';
Height = Vars(:,4*N+1);

%%% Un-normalized depth coordinates
Depth = p.h_0*repmat(Height', N, 1).*repmat(z_h, 1, numel(Time));

zeta830final = interp1(Phi(:,end),Depth/p.h_0,1-830/p.rho_i);  % normalized



%%% Normalized ice velocity
W = nan(N, numel(Time));
Mass = nan(numel(Time),1);
for i = 1:numel(Time)
    S_int = Height(i,1).*(1 - Phi(:,i));
    Sigma = cumtrapz(z_h, S_int);
    V_int = -(Height(i,1)/Ar).*Sigma.^n.*Phi(:,i).^m...
        .*exp(lambda_c*Temp(:,i))./GrainSize(:,i);
    W(:,i) = cumtrapz(z_h, V_int) + w_s;
    
%%% Mass in the column.
    Mass(i) = trapz(Depth(:,i), 1 - Phi(:,i));
end

%%% Create output stucture.
Out = struct('Time', Time, 'Rho', Rho, 'Temp', Temp, 'GrainSize', GrainSize,...
    'Sigma', Sigma, 'Phi', Phi, 'Age', Age, 'Depth', Depth, 'W', W, 'p', p,...
    'H', Height, 'Mass', Mass, 'zeta830final',zeta830final);


function dVarsdt = GovEq(t, Vars)
    phi = Vars(1:N, :);
    r2 = Vars(N+1:2*N, :);
    A = Vars(2*N+1:3*N, :);
    T = Vars(3*N+1:4*N, :);    
    H = Vars(4*N+1, :);
        
%%% Compute gradients.
    dphidz = D1*phi;
    dr2dz = D1*r2;
    dAdz = D1*A;

%%% Compute stress.
    s_int = H*(1 - phi);
    sigma = cumtrapz(z_h, s_int);
    
%%% Compute the ice velocity.
    v_int = -(H/Ar).*sigma.^n.*phi.^m.*exp(lambda_c*T)./r2;
    w = cumtrapz(z_h, v_int) + 1*w_s;

%%% Column height.
    dHdt = w(end) - beta/(1 - phi(end));   
    
%%% Change in porosity.
    dphidt = (1/H)*(D1*((1 - phi).*w) + dHdt*z_h.*dphidz);

%%% Change in square of the grain size.
        if sim_r   % (only if sim_r ==1) 
            dr2dt = (1/H)*(dHdt*z_h - w).*dr2dz + (1 - delta*r2).*exp(lambda_g*T);
        else
            dr2dt = zeros(N,1);
        end


%%% Change in temperature.
    dTdt = zeros(N,1);
    
%%% Change in age.
    dAdt = 1 + (1/H)*(dHdt*z_h - w).*dAdz;
%    dAdt = zeros(N,1);
    
%%% Upper surface boundary conditions.
    dphidt(1) = 0;
    dr2dt(1) = 0;
    dAdt(1) = 0;

%%% Collected ODE vector.    
    dVarsdt = [dphidt; dr2dt; dAdt; dTdt; dHdt];
end

function [position, isterminal, direction] = SteadyState(t, Vars)
    dVdt = GovEq(t, Vars);
%     V_max = max(abs(dVdt(1:1*N)));

    dphidz = D1*Vars(1:N);
    dpdt = dVdt(1:N) - (1/Vars(end))*(dVdt(end)*z_h.*dphidz);
    V_max = max(abs(dpdt));

    detect = V_max < 1e-5;
    position = double(detect);
    isterminal = 1;
    direction = 0;
end

end