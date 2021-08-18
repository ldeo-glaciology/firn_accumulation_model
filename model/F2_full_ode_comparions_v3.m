%% Plot figure 2 for the manuscript 

% Script for comparing the results of the full model to the
% steady state ODE model using different values od the gridspacing dz 

% It plots steady state results from the full model and the ODE in one
% panel, the mismatch between them in another, and the mean and max
% mismatch in an inset. 

% This version of the code uses the setup scripts from Rob Skarbek's
% method-of-lines code.


%% 1. setup axes
clear
figure(2)
clf
tiledlayout(1,2,'TileSpacing','compact')
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% 2. define dz

dz = 0.001:0.001:0.02; 
dz_ref = 0.01;

%% 3. loop over dz
for rr = 1:length(dz)
    
    %%% 3.1 Full model
    %%%% 3.1.1 run full model
    p = FirnSetup3('dz',dz(rr));
    tic;
    Out = FCM9(p);
    RunTime = toc
    fullModelnoT = [Out.Phi(:,end)'; -Out.Sigma(:,end)'; -Out.W(:,end)'; Out.GrainSize(:,end)'; Out.Age(:,end)'];
    Out.z = flip(p.z_h);   
    
    %%%% 3.1.2 plot full model
    if abs(dz(rr) - dz_ref) < 1e-10    % only plot in one iteration
        Figure2fullModelResults = Out;
        save Figure2fullModelResults_v2 Figure2fullModelResults p
        ax1 = nexttile;
        h1 = plot(Out.z,fullModelnoT,'k');              % plot all 6 variables
        text(ax1,-0.13,0.98,'a','units','normalized','FontSize',20)
        xlabel('nondimensional height, $z$','FontSize',12); ylabel('nondimensional variable','FontSize',12)
    end
    
    %%% 3.2 ODE model
    %%%% 3.2.1 run ODE model
    Ar = p.ArthenNumber;
    delta = p.delta;
    phi_s = p.phi_s;
    r_s = p.r2_s_dim/p.r2_0; % grain size at the surface (non-dimensional)
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [zeta,ODEModel] = ode45(@(x,y) fullODEmodel(x,y,Ar,delta),p.z_h,[p.phi_s 0 -p.beta/(1-p.phi_s) r_s 0],options);
    ODEModel = (ODEModel)';    % transpose so that the ODE model results can be compared to the full-model results
    
    
    %%%% 3.2.2 plot ODE model on top of the full model
    if abs(dz(rr) - dz_ref) < 1e-10    % only plot in one iteration
        hold on
        set(gca,'ColorOrderIndex',1)
        h2 = plot(Out.z,ODEModel,'Linewidth',2);
        lgd = legend(h2,'$\phi$','$\sigma$','$w$','$r^2$','$A$','location','southwest','FontSize',12);
    end
    
    %% 3.3 quantitatively compare the results
    
    meanDiff = mean(abs(ODEModel - fullModelnoT),'all');  % the mean absolute difference between the two simulations
    meanDiffPercent = mean(abs((ODEModel - fullModelnoT)./(fullModelnoT)),'all','omitnan')*100;   %  the mean difference as a percentage
    
    maxDiff = max(abs((ODEModel - fullModelnoT)),[],'all');   %  the max absolute difference between the two simulations
    maxDiffPercent = max(abs((ODEModel - fullModelnoT)./(fullModelnoT)),[],'all')*100;   %  the max difference as a percentage
    
    % omit the areas where phi is very small for an alternative calculation of percentage
    I = find(ODEModel(:,1)>1e-3);
    maxDiffPercentLargePhi = max(abs((ODEModel(I,:) - fullModelnoT(I,:))./(fullModelnoT(I,:))),[],'all')*100;   %  the max difference as a percentage
    
    if abs(dz(rr) - dz_ref) < 1e-10    % only in one iteration
        format = 'The mean absolute difference between the two solutions across the five variables is %4.1e and the maximum difference is %4.1e. Mean and maximum differences as a percentage of the full-model values are %3.2f and %3.2f, respectively.';
        sentenceAboutModelMismatch = sprintf(format,meanDiff, maxDiff,meanDiffPercent,maxDiffPercent)      
%         format = 'The model reached steady state in  %4.2f (nondimensional) (%2.0f years)';
%         sentenceAboutReachingSteadyState = sprintf(format,temp.time2ss,temp.time2ss*temp.scales.t0/(365*24*60*60))  
    end
    meanDiffOUT(rr) = meanDiff;
    maxDiffOUT(rr) = maxDiff;
    
    %% 3.4 plot differences in a 2nd panel
 
    if abs(dz(rr) - dz_ref) < 1e-10     % only plot in one iteration
        ax2 = nexttile;
        diffForPlotting = ODEModel-fullModelnoT;
        h3 = plot(Out.z,diffForPlotting,'Linewidth',2);      
        legend(h3,'$\phi$','$\sigma$','$w$','$r^2$','$A$','interpreter','latex','location','southwest','FontSize',12)
        text(ax2,-0.13,0.98,'b','units','normalized','FontSize',20)
        xlabel('nondimensional height, $z$','FontSize',12); ylabel('mismatch between ODE and full-model results','FontSize',12)
        box on
%         title('ODE results - full-model results')
    end
end

%% 4. plot mistmatch as a fn. of dz in an inset
figure(2)
set(gcf,'pos',[ 255   392         883         346])
axes('Position',[0.5968    0.6792    0.1574    0.2027])
h4 = plot(dz,meanDiffOUT,dz,maxDiffOUT);
set(h4,'LineWidth',2) 
set(h4(2),'LineStyle',':') 
xlabel('$\Delta z$','FontSize',12); ylabel('mismatch','FontSize',12)
legend(h4,'mean','maximum','location','northwest')


%% 5. finish figure
h1(1).LineStyle = '-';  h1(1).LineWidth = 3;
h1(2).LineStyle = '--'; h1(2).LineWidth = 3;
h1(3).LineStyle = '-.'; h1(3).LineWidth = 1;
h1(4).LineStyle = '--'; h1(4).LineWidth = 1;
h1(5).LineStyle = '-';  h1(4).LineWidth = 1;

h2(1).LineStyle = '-';  h2(1).LineWidth = 3;
h2(2).LineStyle = '--'; h2(2).LineWidth = 3;
h2(3).LineStyle = '-.'; h2(3).LineWidth = 1;
h2(4).LineStyle = '--'; h2(4).LineWidth = 1;
h2(5).LineStyle = '-';  h2(4).LineWidth = 1;

h3(1).LineStyle = '-';  h3(1).LineWidth = 3;
h3(2).LineStyle = '--'; h3(2).LineWidth = 3;
h3(3).LineStyle = '-.'; h3(3).LineWidth = 1;
h3(4).LineStyle = '--'; h3(4).LineWidth = 1;
h3(5).LineStyle = '-';  h3(4).LineWidth = 1;

set(gcf,'Color','w')
% set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% 6. print figure
print('-dpng','F2_full_ode_comparisons_v5.png')

%% 7. compute some values quoted in the text

%% 8. some numbers for the text

%% How long did it take to reach a steady state?
%%% (99% of its final values
load Figure2fullModelResults_v2 
Out = Figure2fullModelResults; 
z = Out.z;
z0 = Out.p.z_0;
[dt,~] = meshgrid(diff(Out.Time),1:length(z));

I1 = find(all((abs(diff(Out.Phi,1,2))./dt)<1e-3,1),1,'first');
ssTime = Out.Time(I1);
format = 'Porosity approaches a steady state (partial phi / partial t < 10^{-3}) after a nondimensional time of %4.3f (%2.1f years).';
sentenceAboutSteadyState = sprintf(format,ssTime,ssTime*Out.p.t_0/Out.p.spy)  


%% detect the inflection point in the steady state phi profile in the full model. 
% figure (7)
phi = Out.Phi(:,end);
plot(Out.z,phi)

z_mid = z(2:end-1);
I2 = find(z_mid>0.5);
d2phidz2 = diff(phi,2);
inflectionPoint = interp1(d2phidz2(I2),z_mid(I2),0);
% plot(d2phidz2(I2),z_mid(I2),'.-',0,inflectionPoint,'*')
% xlim([-0.00001 0.00001])
format = 'There is an inflection point in phi at z = %4.3f (nondimensional) (%2.1f m below the surface).';
sentenceAboutInflectionPoint = sprintf(format,inflectionPoint,z0-inflectionPoint*z0)  

