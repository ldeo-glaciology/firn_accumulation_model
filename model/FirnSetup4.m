function p = FirnSetup4(inputs)

%% 1. parse input arguments
arguments
    inputs.b0_mpy (1,1) {mustBeNumeric} = 0.1
    inputs.beta (1,1) {mustBeNumeric} = 1
    inputs.beta_multiplier (1,:) {mustBeNumeric} = 1
    inputs.T_s_dim (1,1) {mustBeNumeric} = 253.15
    inputs.plotting (1,1) {mustBeMember(inputs.plotting,[0,1])} = 1
    inputs.saving_xt (1,1) {mustBeMember(inputs.saving_xt,[0,1])} = 1
    inputs.plotting_period (1,1) {mustBeNumeric} = 3000
    inputs.sampling_period (1,1) {mustBeNumeric} = 20
    inputs.dz (1,1) {mustBeNumeric} = 0.01
    inputs.t_total (1,1) {mustBeNumeric} = 20
    inputs.save_dir (1,:) string = 'results_scratch'
    inputs.save (1,1) {mustBeMember(inputs.save,[0,1])} = 1
    inputs.sim_T (1,1) {mustBeMember(inputs.sim_T,[0,1])} = true
    inputs.sim_r (1,1) {mustBeMember(inputs.sim_r,[0,1])} = true
    inputs.PauseGrainEvolution_t (1,1) {mustBeNumeric} = NaN
    inputs.z0 (1,1) {mustBeNumeric} = 100
    inputs.r_s_dim (1,1) {mustBeNumeric} = 2.5000e-07   % grain size at the surface (0.5 mm)^2 
    inputs.phi_s (1,1) {mustBeNumeric} = 0.5   %
    inputs.n (1,1) {mustBeNumeric} = 1
end

%% 2. Define the dimensional parameters of the system.
%%% N D1 beta xi_p delta lambda_g lambda_c m n Pe Ar
%%% These are the dimensional parameters in the problem.
b0_mpy = inputs.b0_mpy;       % ice equivalent accumulation rate [m / yr]
T_s_dim = inputs.T_s_dim;     % upper surface temperature [K]
z_0 = inputs.z0;              % initial column height [m]
dz_dim = inputs.dz*z_0;       % dimensional numerical grid spacing [m]
r2_s_dim = inputs.r_s_dim;    % upper surface grain size [m^2] (0.5 mm)^2 
r2_f = 0.01^2;                % maximum grain size [m^2] (1 cm)^2 
phi_s = inputs.phi_s;         % upper surface porosity
c_i = 2009;                   % heat capacity [J / (kg K)]
E_c = 60e3;                   % compaction activation energy [J]
E_g = 42e3;                   % grain growth activation energy [J]
k_c = 9.2e-9;                 % flow parameter [kg?1 m3 s] from Arthern kc
k_g = 1.3e-7/r2_f;            % grain growth rate constant [m^2 / s]
m = 1;                        % porosity exponent
n = inputs.n;                 % stress exponent
kappa_0 = 2.1;                % thermal conductivity [W / (m K)], % arthern and wingham 1998 pg. 18
G = 0.05;                     % geothermal heat flux  [W/m^2]

%% 3. Constants.
g = 9.81;                           % acceleration due to gravity [m s^-2]
spy = 24*365*3600;                  % seconds per year
R = 8.31;                           % ideal gas constant
rho_i = 918;                        % ice density [kg / m^3]

%% 4. Scales and parameters
%%% Scaling parameters.
b_0 = b0_mpy/spy;
h_0 = z_0;
r2_0 = (h_0*k_g*r2_f/b_0)*exp(-E_g/(R*T_s_dim));
t_0 = h_0/b_0;
T_0 = G*z_0/kappa_0;
sigma_0 = g*rho_i*h_0;
w_0 = b_0;                % scale of the vertical velocity is the accumulation rate

%%% Non-dimensional parameters.
lambda_c = E_c*T_0/(R*T_s_dim^2);
lambda_g = E_g*T_0/(R*T_s_dim^2);
gamma = (sigma_0/4)^(1-n);   % factor to account for non-linear rheology, under the default value of n = 1, gamma = 1 and it has no effect on the value of Ar
Ar = r2_0/(k_c*t_0*sigma_0^n*exp(-E_c/(R*T_s_dim)))/gamma; % Arthern number
Fl = h_0*G/(kappa_0*T_0);
Pe = rho_i*c_i*b_0*h_0/kappa_0;
beta = inputs.beta;
delta = r2_0/r2_f;


%% 5. Set up real space grid in height coordinates.
z_init = (z_0: -dz_dim: 0)';
N = numel(z_init);

%%% Normalized depth coordinates.
z_h = flip(z_init)/z_0;

%% 6. Initial conditions
%%% Initial porosity, temperature, and grain size.
phi_init = (z_init/z_0).*phi_s; % THIS GIVES phi = 0 AAT BOTTOM.
%phi_init = exp(-(flip(z_init)/z_0)).*phi_s;
%phi_init = phi_s*ones(numel(z_init),1);

%%% Compacted grid. Need to flip z_init to be in depth coordinates.
xi_init = flip(z_init) - cumtrapz(flip(z_init), phi_init);

w_s = beta/(1 - phi_s);

%%% Dimensional initial temperature and grain size squared.
T_init = T_0*(1 - z_init/z_0) + T_s_dim;
if inputs.sim_r
    r2_init = r2_s_dim + r2_0.*(1 - z_init/z_0);
else
    r2_init = r2_s_dim + 0*z_init;
end

%%% Dimensionless initial temperature and grain size squared.
T_hat_init = zeros(size(z_init));
%r2_hat_init = 0.01 + (1 - z_init/z_0);
r2_hat_init = r2_init/r2_0;


%%% xi_p is the normalized compacted coordinate.
xi_p = xi_init/xi_init(end);
dxi_p = xi_p(2) - xi_p(1);

%%% Dimensionless initial stress.
sigma_hat_init = cumtrapz(z_h, 1 - phi_init);
%sigma_hat_init = cumsimps(z_h, 1 - phi_init);

%%% Dimensionless firn age.
A_s = 0;
A_hat_init = linspace(A_s, 1, numel(z_init))';

v_int = -(1/Ar).*sigma_hat_init.^n.*phi_init.^m.*exp(lambda_c*T_hat_init)./r2_hat_init;
w_init = cumtrapz(z_h, v_int) + w_s;
%w(end) = beta/(1 - phi(end));
w = flip(cumtrapz(flip(z_h), flip(v_int)) + beta/(1 - phi_init(end)));

%%% Initial conditions vector.
Vars_init = [phi_init; r2_hat_init; sigma_hat_init; A_hat_init; T_hat_init; z_0/h_0];

%% 7. Define gradient operator

%%% Finite difference gradient operator using three point upwind scheme.
%D1 = three_point_upwind_uni_D1(z_h(1), z_h(end), N, 1);
D1 = two_point_upwind_uni_D1(z_h(1), z_h(end), N, 1);
%D1 = three_point_centered_uni_D1(z_h(1), z_h(end), N);
%D2 = three_point_centered_uni_D2(z_h(1), z_h(end), N);

%% 8. Save outputs for use at the solution stage

% Collect all the parameters into a structure 'p'.
p.GridNumber = N;
p.Gradient = D1;
p.spy = spy;
p.beta = beta;
p.xi_p = xi_p;
p.dxi_p = dxi_p;
p.delta = delta;
p.lambda_c = lambda_c;
p.lambda_g = lambda_g;
p.PorosityExponent = m;
p.StressExponent = n;
p.PecletNumber = Pe;
p.FluxNumber = Fl;
p.ArthenNumber = Ar;
p.InitialConditions = Vars_init;
p.rho_i = rho_i;
p.h_0 = h_0;
p.T_0 = T_0;
p.T_s = T_s_dim;
p.t_0 = t_0;
p.b_0 = b_0;
p.r2_0 = r2_0;
p.sigma_0 = sigma_0;
p.r2_s_dim = r2_s_dim;
p.w_0 = w_0;
p.w_s = w_s;
p.z_h = z_h;
p.z_0 = z_0;
p.b0_mpy = b0_mpy;
p.phi_s = phi_s;
p.sim_r = inputs.sim_r;
p.rho_i = rho_i;
p.beta_multiplier = inputs.beta_multiplier;

