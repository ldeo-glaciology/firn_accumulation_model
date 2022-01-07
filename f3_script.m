%% Plots figure 3 for the manuscript.

% a series of experiments with no grain-size evolution. using the full firn
% model and two simplifications to the ODE model.

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


%% 1. setup axes
clear
figure(3)
clf
tiledlayout(3,2,'TileSpacing','none','padding','none')
set(gcf,'pos',[57    84   717   676])
ax1 = nexttile;    ylabel('porosity, $\phi(z)$')
ax2 = nexttile;    ylabel('velocity, $w(z)$')
ax3 = nexttile;    ylabel('porosity, $\phi(z)$')
ax4 = nexttile;    ylabel('velocity, $w(z)$')
ax5 = nexttile;    ylabel('porosity, $\phi(z)$');  xlabel('depth, $z$ ');
ax6 = nexttile;    ylabel('velocity, $w(z)$'); xlabel('depth, $z$ ');

inset1 = axes('Position',[0.2718     0.8121    0.2031    0.1409]);
inset2 = axes('Position',[0.2718     0.505    0.2031    0.1409]); 
inset3 = axes('Position',[0.2718     0.18     0.2031    0.1409]);


h1 = []; h2 = []; h3 = []; h4 = [];h5 = []; h6 = [];


%% 2. define beta vector
% beta
N = 20;
D = 1;
p1 = 0.5; p2 =10; 
beta = p1 + (p2-p1)*linspace(0,1,N).^D;


%% 2. run the full model and two simplifications of the ODE model while varying beta
clear results_1 

rho_i = 918;

for jj = 1:length(beta)
    
    %%% 2.1 full model (r(z,t)=r_s)   
    %%%% 2.1.1 run full model
    p = FirnSetup5('beta',beta(jj),'sim_r',false);
    tic;
    results_1(jj)  = FCM9b(p);
    RunTime = toc
    
    %%%% 2.1.1 pot full model
    hold(ax1,'on')
    finalDepth = p.z_h*results_1(jj).H(end);
    h1 = [h1;  plot(ax1,finalDepth,results_1(jj).Phi(:,end),'k')];
    hold(ax2,'on')
    h2 = [h2; plot(ax2,finalDepth,results_1(jj).W(:,end),'k')];
    
    %%% 2.3 ODE model with \sigma = -zeta
    %%%% 2.3.1 run ODE model 1
    Ar = p.ArthenNumber;
    delta = p.delta;
    phi_s = p.phi_s;
    r_s = p.r2_s_dim/p.r2_0; % grain size at the surface (non-dimensional)   
    n=1; m=1;
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [zeta,y_ODE1] = ode45(@(x,y) odeSig(x,y,n,m,Ar,r_s) ,p.z_h,[phi_s beta(jj)/(1-phi_s)],options);
%     z=flip(zeta);
    
    %%%% 2.3.2 plot ODE model 1
    hold(ax3,'on')
    h3 = [h3; plot(ax3,zeta,y_ODE1(:,1),'k')];
    hold(ax4,'on')
    h4 = [h4; plot(ax4,zeta,y_ODE1(:,2),'k')];
    
    %%%% 2.3.3 compute z830
    zeta830_ODE1(jj) = interp1(y_ODE1(:,1),zeta,1-830/p.rho_i);
    
    %%% 2.4 ODE model 2, with \sigma = -zeta and w = -beta
    %%%% 2.4.1 run ODEmodel 2
    [zeta,y_ODE2] = ode45(@(x,y) odeSigW(x,y,n,m,Ar,r_s,beta(jj)) ,p.z_h,phi_s,options);
    
    %%%% 2.4.2 plot ODE model 2
    hold(ax5,'on')
    h5 = [h5; plot(ax5,zeta,y_ODE2(:,1),'k')];
    hold(ax6,'on')
    h6 = [h6; plot(ax6,zeta,zeta*0+beta(jj),'k')];
    
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
ylim(axs([2 4 6]),[-0.5 20])
xlabel(ax5,'nondimensional depth, $z$')
xlabel(ax6,'nondimensional depth, $z$')

ylabel(inset1,'$z_{830}$')
xlabel(inset1,'$\beta$')
ylabel(inset2,'$z_{830}$')
xlabel(inset2,'$\beta$')
ylabel(inset3,'$z_{830}$')
xlabel(inset3,'$\beta$')

set(ax1,'XTickLabel',[]);
set(ax2,'XTickLabel',[]);
set(ax3,'XTickLabel',[]);
set(ax4,'XTickLabel',[]);

set(findall(gcf,'-property','FontSize'),'FontSize',12)

title(ax1,'Full model, no grain growth','FontSize',15)
title(ax2,'Full model, no grain growth','FontSize',15)
title(ax3,'ODE model (Eqn 25) with $\sigma = -z$','FontSize',15)
title(ax4,'ODE model (Eqn 25) with $\sigma = -z$','FontSize',15)
title(ax5,'ODE model (Eqn 26) with $\sigma = -z$; $w = -\beta$','FontSize',15)
title(ax6,'ODE model (Eqn 26) with $\sigma = -z$; $w = -\beta$','FontSize',15)

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

%% add arrows
annotation(gcf,'textarrow',[0.111576011157601 0.211994421199442],...
    [0.760834319526627 0.801775147928994]);
annotation(gcf,'textarrow',[0.105997210599721 0.179916317991632],...
    [0.417639053254438 0.452662721893491]);
annotation(gcf,'textarrow',[0.099023709902371 0.157601115760112],...
    [0.106988165680473 0.133136094674556]);

annotation(gcf,'textarrow',[0.789 0.789],...
    [0.389532544378698 0.498520710059172]);
annotation(gcf,'textarrow',[0.789 0.789],...
    [0.71201775147929 0.821005917159763]);
annotation(gcf,'textarrow',[0.789 0.789],...
    [0.0700059171597633 0.171597633136095]);

%% 5. print figure
print('-dpng','figures/f3.png')
% 
%% 6. additional values quoted in text
% dimensional grain radius at the surface =
sqrt(p.r2_s_dim)    % 5.0000e-04 m (0.5 mm is quoted in the manuscript).
