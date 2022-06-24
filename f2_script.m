%% Plot figure 2 for the manuscript 

% Script for comparing the results of the full model to the
% steady state ODE model using different values of the gridspacing dz 

% It plots steady state results from the full model and the ODE in one
% panel, the mismatch between them in another, and the mean and max
% mismatch in an inset. 

% This version of the code uses the setup scripts from Rob Skarbek's
% method-of-lines code.

% Modified in response to reviewer comments:
%    -- increased font size on figure axis labels and tick labels

%% 1. setup axes
clear
figure(2)
clf
tiledlayout(1,2,'TileSpacing','compact')
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% 2. define dz

dz = 0.003:0.001:0.02; 

NgridPoints = round(1./dz);

dz = 1./NgridPoints;

dz_ref = dz(8);

%% 3. loop over dz
for rr = 1:length(dz)
    
    %% % 3.1 Full model
    %%%% 3.1.1 run full model
    p = FirnSetup5('dz',dz(rr));

    tic;
    Out = FCM9b(p);

    RunTime = toc
    Out.Age(Out.Age<1e-20) = 0;    % correct for a small numerical error which yields non zero age at the surface (~1e-30)
    fullModelnoT = [Out.Phi(:,end)'; -Out.Sigma(:,end)'; Out.W(:,end)'; Out.GrainSize(:,end)'; Out.Age(:,end)'];

    
    %%%% 3.1.2 plot full model
    if abs(dz(rr) - dz_ref) < 1e-10    % only plot in one iteration
        Figure2fullModelResults = Out;
        save savedResults/f2_results Figure2fullModelResults p
        ax1 = nexttile;
        h1 = plot(p.z_h*Out.H(end),fullModelnoT,'k');              % plot all 6 variables, scaled onto the original z_0 = 100m grid (by multiplying by Out.H(end))
        text(ax1,-0.13,0.98,'a','units','normalized','FontSize',20)
        xlabel(ax1,'nondimensional depth, $z$','FontSize',12); ylabel(ax1,'nondimensional variable','FontSize',12)
    end
    
    %% % 3.2 ODE model
    %%%% 3.2.1 run ODE model
    Ar = p.ArthenNumber;
    delta = p.delta;
    phi_s = p.phi_s;
    r_s = p.r2_s_dim/p.r2_0; % grain size at the surface (non-dimensional)
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [zeta,ODEModel] = ode45(@(x,y) fullODEmodel(x,y,Ar,delta),p.z_h,[p.phi_s 0 p.beta/(1-p.phi_s) r_s 0],options);
    ODEModel = ODEModel';    % transpose so that the ODE model results can be compared to the full-model results
    
    
    %%%% 3.2.2 plot ODE model on top of the full model
    if abs(dz(rr) - dz_ref) < 1e-10    % only plot in one iteration
        hold on
        set(gca,'ColorOrderIndex',1)
        h2 = plot(p.z_h,ODEModel,'Linewidth',2);
        lgd = legend(h2,'$\phi$','$\sigma$','$w$','$r^2$','$A$','location','southwest','FontSize',12);
    end
    
    %% % 3.3 quantitatively compare the results
    %%%% 3.3.1 interp one of the results onto the other (depending on which grid covers a larger dimensional distance).
    % The ODE model grid always covers z0 (100m), so Out.H(end) controls if
    % the full model exceeds this or not. 
    
    if Out.H(end) <= 1       % interp ODE results onto full-model grid
        ODEModel_new = interp1(p.z_h,ODEModel',p.z_h*Out.H(end))';
        fullModelnoT_new = fullModelnoT;   % unchanged
        z_for_comparison = p.z_h*Out.H(end);
    elseif Out.H(end) > 1    % interp full results onto ODE-model grid
        error('code only currently works when H<1')
    end
        
    meanDiff = mean(abs(ODEModel_new - fullModelnoT_new),'all');  % the mean absolute difference between the two simulations
    meanDiffPercent = mean(abs((ODEModel_new - fullModelnoT_new)./(fullModelnoT_new)),'all','omitnan')*100;   %  the mean difference as a percentage
    
    maxDiff = max(abs((ODEModel_new - fullModelnoT_new)),[],'all');   %  the max absolute difference between the two simulations
    maxDiffPercent = max(abs((ODEModel_new - fullModelnoT_new)./(fullModelnoT_new)),[],'all')*100;   %  the max difference as a percentage
    
    % omit the areas where phi is very small for an alternative calculation of percentage mismatch
    I = find(ODEModel_new(1,:)>1e-2);
    maxDiffPercentLargePhi = max(abs((ODEModel_new(:,I) - fullModelnoT_new(:,I))./(fullModelnoT_new(:,I))),[],'all')*100;   %  the max difference as a percentage
    
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
        diffForPlotting = ODEModel_new-fullModelnoT_new;
        h3 = plot(z_for_comparison,diffForPlotting,'Linewidth',2);      
        legend(h3,'$\phi$','$\sigma$','$w$','$r^2$','$A$','interpreter','latex','location','southwest','FontSize',12)
        text(ax2,-0.13,0.98,'b','units','normalized','FontSize',20)
        xlabel(ax2,'nondimensional depth, $z$','FontSize',12); ylabel(ax2,'mismatch between ODE and full-model results','FontSize',12)
        box on
%         title('ODE results - full-model results')
    end
end

%% 4. plot mistmatch as a fn. of dz in an inset
figure(2)
set(gcf,'pos',[ 255         291        1025         447])
axes('Position',[0.7548    0.1879    0.1311    0.1712])
h4 = plot(dz,meanDiffOUT,dz,maxDiffOUT);
set(h4(1),'LineWidth',3,'Color',[0.5 0.5 0.5]) 
set(h4(2),'LineWidth',2,'LineStyle',':') 
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


% increase size of tick labels
ax1.XAxis.FontSize = 13;
ax1.XLabel.FontSize = 16;
ax1.YAxis.FontSize = 13;
ax1.YLabel.FontSize = 16;
ax2.XAxis.FontSize = 13;
ax2.XLabel.FontSize = 16;
ax2.YAxis.FontSize = 13;
ax2.YLabel.FontSize = 16;


set(gcf,'Color','w')
% set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% 6. print figure
print('-dpng','figures/f2.png')

%% 7. compute some values quoted in the text


%%% How long did it take to reach a steady state?
%%% (99% of its final values)
load savedResults/f2_results 
Out = Figure2fullModelResults; 
% z = Out.z;
% z0 = Out.p.z_0;
[dt,~] = meshgrid(diff(Out.Time),1:size(Out.Depth,1));

I1 = find(all((abs(diff(Out.Phi,1,2))./dt)<1e-3,1),1,'first');
ssTime = Out.Time(I1);
format = 'Porosity approaches a steady state (partial phi / partial t < 10^{-3}) after a nondimensional time of %4.3f (%2.1f years).';
sentenceAboutSteadyState = sprintf(format,ssTime,ssTime*Out.p.t_0/Out.p.spy)  


%%% detect the inflection point in the steady state phi profile in the full model. 
% figure (7)
finalPhi = Out.Phi(:,end);
finalDepth = Out.Depth(:,end)/p.h_0;
% plot(finalDepth,finalPhi)

finalDepth_mid = finalDepth(2:end-1);
I2 = find(finalDepth_mid<0.5);
d2phidz2 = diff(finalPhi,2);
inflectionPoint = interp1(d2phidz2(I2),finalDepth_mid(I2),0);
% hold on
% plot(finalDepth_mid,d2phidz2,'.-',inflectionPoint,0,'*')
% axis([0.1477    0.2793   -0.0012    0.0012])
format = 'There is an inflection point in phi at zeta = %4.3f (nondimensional) (%2.1f m).';
sentenceAboutInflectionPoint = sprintf(format,inflectionPoint,inflectionPoint*p.h_0)  

