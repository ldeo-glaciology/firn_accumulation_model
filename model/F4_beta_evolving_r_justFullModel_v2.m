%% Plots figure 4 for the manuscript.

% A series of experiments with the non-dimensional firn model with grain
% size evolution. It runs the model for multiple combinations of beta and r_s for the full
% contour plot in figure 6a and saves the results then plots three
% sets of results. 

% This version of the code uses the setup scripts from Rob Skarbek's
% method-of-lines code.

%% 1. setup axes
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


clear
figure(4)
clf
tiledlayout(3,2,'TileSpacing','none','padding','none')
set(gcf,'pos',[57    84   717   676])
ax1 = nexttile;  
ax2 = nexttile;  
ax3 = nexttile;  
ax4 = nexttile;  
ax5 = nexttile;   
ax6 = nexttile;    

inset1 = axes('Position',[0.11     0.826   0.2031    0.1409]);
inset2 = axes('Position',[0.11    0.505   0.2031    0.1409]);
inset3 = axes('Position',[0.11     0.18    0.2031   0.1409]);

%% 2. define beta and r_S
% beta
N = 21;
D = 3;
p1 = 0.1; p2 =10; 
beta = p1 + (p2-p1)*linspace(0,1,N).^D;

% r_s
p = FirnSetup3; % run the setup script to get the grain size scale, r2_0
r0 = p.r2_0;
N = 21;    % define a range of the non-dimensional surface grain size to vary over
D = 3;
p1 = 0.001; p2 = 0.1; 
r_s  = p1 + (p2-p1)*linspace(0,1,N).^D;
r_s_dim = r_s*r0;   % change this into a dimensional surface grain size using, r2_0

%% 3. run the full model while thickening
load rerun rerun
if rerun
    for ii  = 1:length(r_s_dim)
        for jj = 1:length(beta)
            p = FirnSetup3('beta',beta(jj),'r_s_dim',r_s_dim(ii));
            results_1(jj,ii)  = FCM9(p);
        end
        ii-length(r_s_dim)
    end
 % save the results for use later
 reduceResultsStruct     % remove all but the final timestep in each of the variables, using this short script. 
 save results_for_Figure4_4 results_1 r_s r0 r_s_dim beta
end

%% 4. reload the results of this loop, to save time when replotting
% (to rerun the results, uncomment the cell above and comment out the load
% command below, then rerun the whole script.) 
if ~rerun
    load results_for_Figure4_4
end

%% 5. plot three sets of results with the same r_s
I = [1 10 21];

% initiate empty arrays for the figure handles. 
h1 = []; h2 = []; h3 = []; h4 = []; h5 = []; h6 = [];


for kk = 1:length(beta)   %  plot three sets of results
    
    % plot top row (lowest r_s)           
    hold(ax1,'on')
    h1 = [h1; plot(ax1,[results_1(kk,I(1)).z],[results_1(kk,I(1)).Phi],'k')];
    hold(ax2,'on')
    h2 = [h2; plot(ax2,[results_1(kk,I(1)).z],-[results_1(kk,I(1)).W],'k')];    
   

    % plot second row (middle r_s)           
    hold(ax3,'on')
    h3 = [h3; plot(ax3,[results_1(kk,I(2)).z],[results_1(kk,I(2)).Phi],'k')];
    hold(ax4,'on')
    h4 = [h4; plot(ax4,[results_1(kk,I(2)).z],-[results_1(kk,I(2)).W],'k')];
    
    % plot bottom row (largest r_s)         
    hold(ax5,'on')
    h5 = [h5; plot(ax5,[results_1(kk,I(3)).z],[results_1(kk,I(3)).Phi],'k')];
    hold(ax6,'on')
    h6 = [h6; plot(ax6,[results_1(kk,I(3)).z],-[results_1(kk,I(3)).W],'k')];
end


%% 6. plot the insets
% inset1
hold(inset1,'on')
for ii  = 1:2:length(r_s_dim)
    plot(inset1,beta,[results_1(:,ii).zeta830final],'.:k')
end
plot(inset1,beta,[results_1(:,I(1)).zeta830final],'.-k','LineWidth',1)

% inset2
hold(inset2,'on')
for ii  = 1:2:length(r_s_dim)
    plot(inset2,beta,[results_1(:,ii).zeta830final],'.:k')
end
plot(inset2,beta,[results_1(:,I(2)).zeta830final],'.-k','LineWidth',1)

% inset3
hold(inset3,'on')
for ii  = 1:2:length(r_s_dim)
    plot(inset3,beta,[results_1(:,ii).zeta830final],'.:k')
end
plot(inset3,beta,[results_1(:,I(3)).zeta830final],'.-k','LineWidth',1)

%% 7. finish plots
axs = [ax1; ax2; ax3; ax4; ax5; ax6];
ylim(axs([1 3 5]),[0 0.5])
ylim(axs([2 4 6]),[-20 0.5])
ylim([inset1; inset2; inset3],[0 1])

ylabel(ax1,'porosity, $\phi(z)$')
ylabel(ax2,'nondimensional velocity, $w(z)$')
ylabel(ax3,'porosity, $\phi(z)$')
ylabel(ax4,'nondimensional velocity, $w(z)$')
ylabel(ax5,'porosity, $\phi(z)$')
ylabel(ax6,'nondimensional velocity, $w(z)$')

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

title(ax1,['Full model, evolving grain-size: $r_s^2$ = ' num2str(r_s(I(1)),'%4.3f')],'FontSize',15)
title(ax2,['Full model, evolving grain-size: $r_s^2$ = ' num2str(r_s(I(1)),'%4.3f')],'FontSize',15)
title(ax3,['Full model, evolving grain-size: $r_s^2$ = ' num2str(r_s(I(2)),'%4.2f')],'FontSize',15)
title(ax4,['Full model, evolving grain-size: $r_s^2$ = ' num2str(r_s(I(2)),'%4.2f')],'FontSize',15)
title(ax5,['Full model, evolving grain-size: $r_s^2$ = ' num2str(r_s(I(3)),'%4.1f')],'FontSize',15)
title(ax6,['Full model, evolving grain-size: $r_s^2$ = ' num2str(r_s(I(3)),'%4.1f')],'FontSize',15)

text(ax1,-0.11,1,'a','units','normalized','FontSize',20)
text(ax2,-0.11,1,'b','units','normalized','FontSize',20)
text(ax3,-0.11,1,'c','units','normalized','FontSize',20)
text(ax4,-0.11,1,'d','units','normalized','FontSize',20)
text(ax5,-0.11,1,'e','units','normalized','FontSize',20)
text(ax6,-0.11,1,'f','units','normalized','FontSize',20)

% reposition insets
inset1.Position(2) = ax1.Position(2)+0.12;
inset2.Position(2) = ax3.Position(2)+0.12;
inset3.Position(2) = ax5.Position(2)+0.12;

set(gcf,'Color','w')

ax1.Box = 1;
ax2.Box = 1;
ax3.Box = 1;
ax4.Box = 1;
ax5.Box = 1;
ax6.Box = 1;

%% 8. print figure
print('-dpng','F4_beta_evolving_rJustFull_v2.png')

%% 9. compute some values quoted in the text

% figure(1234)
% load results_for_Figure4_4  

%%% for the higher r_s
N=19;
% plot(beta,[results_1(:,N).zeta830final],'-r','LineWidth',2)
p = polyfit(beta,[results_1(:,N).zeta830final],1);
% hold on
% xplot = 0:0.1:10;
% plot(xplot,xplot*p(1)+p(2))

disp(['mean accumulation dependence when r_s = ' num2str(r_s(N)) ' is ' num2str(p(1))])
%%%% compute the gradient as a function of beta
G = gradient([results_1(:,N).zeta830final],beta);
disp(['(varying from ' num2str(G(1)) ' at $\beta$ = ' num2str(beta(1)) ])
disp([' to ' num2str(G(end)) ' at $\beta$ = ' num2str(beta(end)) ')'])

% (results_1(2,N).zeta830final - results_1(1,N).zeta830final) /(beta(2)-beta(1))   % just as a check that gradient.m is working correctly, this equals G(1)


%%% for the lower r_s
N=1;
% plot(beta,[results_1(:,N).zeta830final],'-g','LineWidth',2)
p = polyfit(beta,[results_1(:,N).zeta830final],1);
% plot(xplot,xplot*p(1)+p(2))
disp(['mean accumulation dependence when r_s = ' num2str(r_s(N)) ' is ' num2str(p(1))])


%%%% compute the gradient
G = gradient([results_1(:,N).zeta830final],beta);
disp(['(varying from ' num2str(G(1)) ' at $\beta$ = ' num2str(beta(1)) ])
disp([' to ' num2str(G(end)) ' at $\beta$ = ' num2str(beta(end))  ')'])
