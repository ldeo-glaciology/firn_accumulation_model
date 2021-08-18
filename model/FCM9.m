function Out = FCM9(p)

N = p.GridNumber;
D1 = p.Gradient;
z_h = p.z_h;
beta = p.beta;
delta = p.delta;
lambda_c = p.lambda_c;
lambda_g = p.lambda_g;
m = p.PorosityExponent;
n = p.StressExponent;
Pe = p.PecletNumber;
Fl = p.FluxNumber;
Ar = p.ArthenNumber;
w_s = p.w_s;
Vars_init = p.InitialConditions;
sim_r = p.sim_r;
beta_multiplier = p.beta_multiplier;

%%% Solve the system
time = [0 10];

options = odeset('Events', @SteadyState, 'RelTol',1e-8,'AbsTol', 1e-8,'InitialStep',1e-3);
[Time, Vars] = ode45(@GovEq, time, Vars_init, options);

%%% Collect the simulation variables.
Phi = Vars(:, 1:N)';
Rho = p.rho_i*(1 - Phi);
GrainSize = Vars(:, N+1:2*N)';
Sigma = Vars(:, 2*N+1:3*N)';
Age = Vars(:, 3*N+1:4*N)';
Temp = Vars(:, 4*N+1:5*N)';

% scaled height (1 is the surface)
z = flip(p.z_h);

zeta830final = interp1(Phi(:,end),p.z_h,1-830/p.rho_i);
deltaAge = interp1(Phi(:,end),Age(:,end),1-830/p.rho_i);

%%% Normalized ice velocity
W = nan(N, numel(Time));
V_int = -(1/Ar).*Sigma.^n.*Phi.^m.*exp(lambda_c*Temp)./GrainSize;
for i = 1:numel(Time)
    W(:,i) = cumtrapz(z_h, V_int(:,i)) + w_s;
end

%%% Un-normalized depth coordinates
Depth = p.h_0*p.z_h;

%%% Create output stucture.
Out = struct('Time', Time, 'Rho', Rho, 'Temp', Temp, 'GrainSize', GrainSize,...
    'Sigma', Sigma, 'Phi', Phi, 'Age', Age, 'Depth', Depth, 'W', W, 'p', p,...
    'zeta830final',zeta830final,'deltaAge',deltaAge,'z',z);

    function dVarsdt = GovEq(t, Vars)
        phi = Vars(1:N, :);
        r2 = Vars(N+1:2*N, :);
        sigma = Vars(2*N+1:3*N, :);
        A = Vars(3*N+1:4*N, :);
        T = Vars(4*N+1:5*N, :);
        
        %%% Compute spatial gradients.
        DR2_1 = D1*r2;
        DA_1 = D1*A;
        
        %%% Compute the ice velocity.
        v_int = -(1/Ar).*sigma.^n.*phi.^m.*exp(lambda_c*T)./r2;
        w = cumtrapz(z_h, v_int) + w_s;
        
        %%% Change in stress.
        dsigmadt = (1 - phi(1))*w_s - (1 - phi).*w;
        
        %%% Change in porosity.
        dphidt = D1*((1 - phi).*w);
        
        %%% Change in square of the grain size.
        if sim_r   % (only if sim_r ==1)
            dr2dt = -w.*DR2_1 + (1 - delta*r2).*exp(lambda_g*T);
        else
            dr2dt = zeros(N,1);
        end
        
        %%% Change in temperature.
        dTdt = zeros(N,1);
        
        %%% Change in age.
        dAdt = 1 - w.*DA_1;
        
        %%% Upper surface boundary conditions.
        dphidt(1) = 0;
        dr2dt(1) = 0;
        dAdt(1) = 0;
        
        %%% Collected ODE vector.
        dVarsdt = [dphidt; dr2dt; dsigmadt; dAdt; dTdt];
    end

    function [position, isterminal, direction] = SteadyState(t, Vars)
        dVdt = GovEq(t, Vars);
        V_max = max(abs(dVdt(1:4*N)));
        detect = V_max < 1e-2;
        position = double(detect);
        isterminal = 1;
        direction = 0;
    end

end