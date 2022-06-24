%% Plots figure 6 for the manuscript  

% Produces contour plots of firn thickness, zeta_830, as a function of
% accumulation rate, beta, and surface grain size, r_s.

% One panel is from the full model. This script loads the results of
% running the model many times with many combinations of beta and r_s
% (produced using beta_evolving_r_justFullModel_v2.m)

% The other two panels are from two different simplified ODE models, which 
% both have evolving grain size: (odeSig_r and odeSigW_r)

% Modified in response to reviewer comments:
%    -- increased font size on figure axis labels and tick labels

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
% F4_beta_evolving_r_justFullModel_v3.m with rerun = 1
load savedResults/f4_results

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
        [zeta,y_ODE1] = ode45(@(x,y) odeSig_r(x,y,n,m,Ar) ,results_1(jj,ii).p.z_h,[phi_s beta(jj)/(1-phi_s) ,r_s(ii)],options);
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
% add panel labels
text(ax1,-0.15,0.98,'a','units','normalized','FontSize',20)
text(ax2,-0.15,0.98,'b','units','normalized','FontSize',20)
text(ax3,-0.15,0.98,'c','units','normalized','FontSize',20)

% add axis labels
ylabel(ax1,'surface grain size, $r^2_s$')
xlabel(ax1,'accumulation rate, $\beta$')
ylabel(ax2,'surface grain size, $r^2_s$')
xlabel(ax2,'accumulation rate, $\beta$')
ylabel(ax3,'surface grain size, $r^2_s$')
xlabel(ax3,'accumulation rate, $\beta$')

% add titles
title(ax1,'full model','FontSize',18)
title(ax2,'ODE model (Eqn 26) with $\sigma = -z$','FontSize',18)
title(ax3,'ODE model (Eqn 27) with $\sigma = -z; w = \beta$','FontSize',18)

% make the color scale the same across panels and figures
caxis(ax1,[0 1])
caxis(ax2,[0 1])
caxis(ax3,[0 1])

% change size of tick labels
ax1 = axis_font_sizes(ax1,13,17);
ax2 = axis_font_sizes(ax2,13,17);
ax3 = axis_font_sizes(ax3,13,17);

%% 6. print figure
print('-dpng','figures/f6.png')
