%% Scripts that plot the figures for the manuscript

%%% to re-run all the simulations set rerun = 1
% if rerun = 0,   in some cases saved results will be reloaded to save
% time. 

rerun = 0
save rerun rerun
close all
%% figure 1 and the scale and parameters table:

F1_ParametersTableAndFigure_v2    % updated tidied
 
%% figure 2 (comparison between full model and ODE model results)

F2_full_ode_comparions_v3    % updated tidies

%% figure 3 (dependence on accumulation rate with no grain size evolution for the full model and ODE model)

F3_beta_NonEvolving_r_v2      % updated  tidied

%% figure 4 (dependence on accumulation rate with grain size evolution for the full model)

F4_beta_evolving_r_justFullModel_v2   % updated  tidied

%% figure 5 (dependence on accumulation rate with grain size evolution for the full model 2: gradient plot)

F5_zeta830_grad_plot_v2   % updated  tidied

%% figure 6    (contour plots of firn thickness for full model and two ODE models)

F6_contours_all_models_v2  % updated  tidied

%% figure 7 (contour plots for three nonlinear rheologys, n ~= 1)

F7_nonlinear_rheology_v2   % updated   tidied

%% figure 8 (contour plots of thickening and thinning experiments). 

% F8_thin_thick_plotting

