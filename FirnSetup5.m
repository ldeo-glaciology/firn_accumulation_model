function p = FirnSetup5(inputs)

%{

This function generates the inputs for FCM9b.m, which solves the model
equations. 

Use FirnSetup5.m as follows:

    1. To use default values

    >> p = FirnSetup5

    in this case all the default valuesdefined in the arguments block below
    will be used. 

    2. To prescribe one or more inputs use name-value pairs 

    >> p = FirnSetup5('dz',0.01,'beta',0.3)


Rob Skarbek, Jonny Kingslake, Lamont-Doherty Earth Observatory, 2021-22

%}
%% 1. Parse input arguments
arguments
    inputs.b0_mpy (1,1) {mustBeNumeric} = 0.1
    inputs.beta (1,1) {mustBeNumeric} = 1
    inputs.nu (1,:) {mustBeNumeric} = 1
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
    inputs.simDuration (1,1) {mustBeNumeric} = 10
    inputs.scaleDuration (1,1) {mustBeNumeric} = 0
end

%% 2. Define the dimensional parameters of the system.
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
nu = inputs.nu;      % accumulation multiplier

%% 5. compute vertical velocity upper biundary conditon
w_s = nu*beta/(1 - phi_s); 

%% 6. Set up real space grid in height coordinates.
z_init = (z_0: -dz_dim: 0)';
N = numel(z_init);

%%% Normalized depth coordinates.
z_h = flip(z_init)/z_0;

%% 7. Initial conditions
%%% Initial porosity
phi_init = (1-z_h).*phi_s;

%%% Dimensional  grain size squared.
if inputs.sim_r
    r2_hat_init = r2_s_dim/r2_0 + z_h;
else
    r2_hat_init = r2_s_dim/r2_0 + 0*z_h;
end

%%% Dimensionless temperature (mostly not used)
T_hat_init = zeros(size(z_init));

%%% Dimensionless firn age.
A_hat_init = z_h;

%%% Initial conditions vector.
Vars_init = [phi_init; r2_hat_init; A_hat_init; T_hat_init; z_0/h_0];

%% 8. Define gradient operator
%%% Finite difference gradient operator using three point upwind scheme.
D1 = two_point_upwind_uni_D1(z_h(1), z_h(end), N, 1);

%% 9. Save outputs for use at the solution stage

% Collect all the parameters into a structure 'p'.
p.GridNumber = N;
p.Gradient = D1;
p.spy = spy;
p.beta = beta;
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
p.nu = nu;
p.simDuration = inputs.simDuration;
p.scaleDuration = inputs.scaleDuration;

