%% Plots figure 3 for the manuscript.

% a series of experiments with no grain-size evolution. using the full firn
% model and two simplifications to the ODE model.

% This version of the code uses the setup scripts from Rob Skarbek's
% method-of-lines code.


warning('off', 'MATLAB:MKDIR:DirectoryExists');
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


%% 1. setup axes
clear
figure(3)
clf
tiledlayout(3,2,'TileSpacing','none','padding','none')
set(gcf,'pos',[57    84   717   676])
ax1 = nexttile;    ylabel('porosity, $\phi(\zeta)$')
ax2 = nexttile;    ylabel('velocity, $w(\zeta)$')
ax3 = nexttile;    ylabel('porosity, $\phi(\zeta)$')
ax4 = nexttile;    ylabel('velocity, $w(\zeta)$')
ax5 = nexttile;    ylabel('porosity, $\phi(\zeta)$');  xlabel('depth, $\zeta$ ');
ax6 = nexttile;    ylabel('velocity, $w(\zeta)$'); xlabel('depth, $\zeta$ ');

inset1 = axes('Position',[0.11     0.826   0.2031    0.1409]);
inset2 = axes('Position',[0.11    0.505   0.2031    0.1409]);
inset3 = axes('Position',[0.11     0.18    0.2031   0.1409]);


h1 = []; h2 = []; h3 = []; h4 = [];h5 = []; h6 = [];


%% 2. define beta vector
% beta
N = 20;
D = 3;
p1 = 0.5; p2 =10; 
beta = p1 + (p2-p1)*linspace(0,1,N).^D;


%% 2. run the full model and two simplifications of the ODE model while varying beta
clear results_1 

rho_i = 918;

for jj = 1:length(beta)
    
    %%% 2.1 full model (r(z,t)=r_s)   
    %%%% 2.1.1 run full model
    p = FirnSetup3('beta',beta(jj),'sim_r',false);
    tic;
    results_1(jj)  = FCM9(p);
    RunTime = toc;
    
    %%%% 2.1.1 pot full model
    hold(ax1,'on')
    h1 = [h1;  plot(ax1,results_1(jj).z,results_1(jj).Phi(:,end),'k')];
    hold(ax2,'on')
    h2 = [h2; plot(ax2,results_1(jj).z,-results_1(jj).W(:,end),'k')];
    
    %%% 2.3 ODE model with \sigma = -zeta
    %%%% 2.3.1 run ODE model 1
    Ar = p.ArthenNumber;
    delta = p.delta;
    phi_s = p.phi_s;
    r_s = p.r2_s_dim/p.r2_0; % grain size at the surface (non-dimensional)   
    n=1; m=1;
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [zeta,y_ODE1] = ode45(@(x,y) odeSig(x,y,n,m,Ar,r_s) ,p.z_h,[phi_s -beta(jj)/(1-phi_s)],options);
    z=flip(zeta);
    
    %%%% 2.3.2 plot ODE model 1
    hold(ax3,'on')
    h3 = [h3; plot(ax3,z,y_ODE1(:,1),'k')];
    hold(ax4,'on')
    h4 = [h4; plot(ax4,z,y_ODE1(:,2),'k')];
    
    %%%% 2.3.3 compute z830
    zeta830_ODE1(jj) = interp1(y_ODE1(:,1),zeta,1-830/p.rho_i);
    
    %%% 2.4 ODE model 2, with \sigma = -zeta and w = -beta
    %%%% 2.4.1 run ODEmodel 2
    [zeta,y_ODE2] = ode45(@(x,y) odeSigW(x,y,n,m,Ar,r_s,beta(jj)) ,p.z_h,phi_s,options);
    
    %%%% 2.4.2 plot ODE model 2
    hold(ax5,'on')
    h5 = [h5; plot(ax5,z,y_ODE2(:,1),'k')];
    hold(ax6,'on')
    h6 = [h6; plot(ax6,z,z*0-beta(jj),'k')];
    
    %%%% 2.4.3 compute z830
    zeta830_ODE2(jj) = interp1(y_ODE2(:,1),zeta,1-830/p.rho_i);
    
end

%% 3. plot the insets

% inset 1
hold(inset1,'on')
plot(inset1,beta,[results_1.zeta830final],'.-k','LineWidth',1)
plot(inset1,beta,zeta830_ODE1,'.:k')
plot(inset1,beta,zeta830_ODE2,'.:k')

% inset 2
hold(inset2,'on')
plot(inset2,beta,[results_1.zeta830final],'.:k')
plot(inset2,beta,zeta830_ODE1,'.-k','LineWidth',1)
plot(inset2,beta,zeta830_ODE2,'.:k')

% inset3
hold(inset3,'on')
plot(inset3,beta,[results_1.zeta830final],'.:k')
plot(inset3,beta,zeta830_ODE1,'.:k')
plot(inset3,beta,zeta830_ODE2,'.-k','LineWidth',1)

%% 4. finish plots
axs = [ax1; ax2; ax3; ax4; ax5; ax6];
ylim(axs([1 3 5]),[0 0.5])
ylim(axs([2 4 6]),[-20 0.5])
ylabel(ax1,'porosity, $\phi(z)$')
ylabel(ax3,'porosity, $\phi(z)$')
ylabel(ax5,'porosity, $\phi(z)$')
xlabel(ax5,'nondimensional height, $z$ ')
xlabel(ax6,'nondimensional height, $z$ ')

ylabel(inset1,'$\zeta_{830}$')
xlabel(inset1,'$\beta$')
ylabel(inset2,'$\zeta_{830}$')
xlabel(inset2,'$\beta$')
ylabel(inset3,'$\zeta_{830}$')
xlabel(inset3,'$\beta$')

set(ax1,'XTickLabel',[]);
set(ax2,'XTickLabel',[]);
set(ax3,'XTickLabel',[]);
set(ax4,'XTickLabel',[]);

set(findall(gcf,'-property','FontSize'),'FontSize',12)

title(ax1,'Full model, no grain growth','FontSize',15)
title(ax2,'Full model, no grain growth','FontSize',15)
title(ax3,'ODE model (Eqn 25) with $\sigma = -\zeta$','FontSize',15)
title(ax4,'ODE model (Eqn 25) with $\sigma = -\zeta$','FontSize',15)
title(ax5,'ODE model (Eqn 26) with $\sigma = -\zeta$; $w = -\beta$','FontSize',15)
title(ax6,'ODE model (Eqn 26) with $\sigma = -\zeta$; $w = -\beta$','FontSize',15)

text(ax1,-0.11,0.98,'a','units','normalized','FontSize',20)
text(ax2,-0.11,0.98,'b','units','normalized','FontSize',20)
text(ax3,-0.11,0.98,'c','units','normalized','FontSize',20)
text(ax4,-0.11,0.98,'d','units','normalized','FontSize',20)
text(ax5,-0.11,0.98,'e','units','normalized','FontSize',20)
text(ax6,-0.11,0.98,'f','units','normalized','FontSize',20)

% reposition insets
% set(inset2,'Position',(ax2.Position(2)+ax2.Position(4))-0.01)
inset1.Position(2) = ax1.Position(2)+0.12;
inset2.Position(2) = ax3.Position(2)+0.12;
inset3.Position(2) = ax5.Position(2)+0.12;

set(inset1,'Color','none')

ax1.Box = 1;
ax2.Box = 1;
ax3.Box = 1;
ax4.Box = 1;
ax5.Box = 1;
ax6.Box = 1;

%% 5. print figure
print('-dpng','F3_beta_NonEvolving_r_v2.png')

%% 6. additional values quoted in text
% dimensional grain radius at the surface =
sqrt(p.r2_s_dim)    % 5.0000e-04 m (0.5 mm is quoted in the manuscript).
