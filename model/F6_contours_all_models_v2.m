%% Plots figure 6 for the manuscript  

% Produces contour plots of firn thickness, zeta_830, as a function of
% accumulation rate, beta, and surface grain size, r_s.

% One panel is from the full model. This script loads the results of
% running the model many times with many combinations of beta and r_s
% (produced using beta_evolving_r_justFullModel_v2.m)

% The other two panels are from two different simplified ODE models, which 
% both have evolving grain size: (odeSig_r and odeSigW_r)

% This version of the code uses the setup scripts from Rob Skarbek's
% method-of-lines code.

%% 1. setup axes
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

clear
figure(6)
clf
tiledlayout(1,3,'TileSpacing','none','padding','none')
set(gcf,'pos',[ 57         436        1069         324])
ax1 = nexttile;
ax2 = nexttile;
ax3 = nexttile;

%% 2. load the parameter-search results
% load the results of running the full model many times for many values
% of r_s and beta - to produce these results run
% F4_beta_evolving_r_justFullModel_v2.m with rerun = 1
load results_for_Figure4_4 

phi_s = results_1(1).p.phi_s;
rho_i = results_1(1).p.rho_i;

%% 3. loop over every value of beta and r_s
tic
for  ii = 1:length(r_s)
    for jj = 1:length(beta)
        % extract the full model results
        zeta_830_full(jj,ii) = results_1(jj,ii).zeta830final;
        
        % 3-equation model (eqn 27)
        Ar = results_1(1).p.ArthenNumber;         % it doesnt matter which one we choose because they are all the same.
        n=1; m=1;
        options = odeset('RelTol',1e-10,'AbsTol',1e-10);
        [zeta,y_ODE1] = ode45(@(x,y) odeSig_r(x,y,n,m,Ar) ,results_1(jj,ii).p.z_h,[phi_s -beta(jj)/(1-phi_s) ,r_s(ii)],options);
        zeta830_ODE1(jj,ii) = interp1(y_ODE1(:,1),zeta,1-830/rho_i);
     
        % 1-equation model (eqn 28)
        [zeta,y_ODE2] = ode45(@(x,y) odeSigW_r(x,y,n,m,r_s(ii),Ar,beta(jj)) ,results_1(jj,ii).p.z_h,phi_s,options);
        zeta830_ODE2(jj,ii) = interp1(y_ODE2(:,1),zeta,1-830/rho_i);

    end
end
toc

%% 4. plot three contour maps: full model and two ode models
contourf(ax1,beta,r_s,zeta_830_full','ShowText','on')
contourf(ax2,beta,r_s,zeta830_ODE1','ShowText','on')
contourf(ax3,beta,r_s,zeta830_ODE2','ShowText','on')

%% 5. finish figure
text(ax1,-0.11,0.98,'a','units','normalized','FontSize',20)
text(ax2,-0.11,0.98,'b','units','normalized','FontSize',20)
text(ax3,-0.11,0.98,'c','units','normalized','FontSize',20)

ylabel(ax1,'surface grain size, $r^2_s$','FontSize',15)
xlabel(ax1,'accumulation rate, $\beta$','FontSize',15)
ylabel(ax2,'surface grain size, $r^2_s$','FontSize',15)
xlabel(ax2,'accumulation rate, $\beta$','FontSize',15)
ylabel(ax3,'surface grain size, $r^2_s$','FontSize',15)
xlabel(ax3,'accumulation rate, $\beta$','FontSize',15)

title(ax1,'full model','FontSize',15)
title(ax2,'ODE model (Eqn 27) with $\sigma = -\zeta$','FontSize',15)
title(ax3,'ODE model (Eqn 28) with $\sigma = -\zeta; w = -\beta$','FontSize',15)

caxis(ax2,[0 1])
caxis(ax1,[0 1])
caxis(ax3,[0 1])

%% 6. print figure
print('-dpng','F6_coutours_all_models_v3.png')
