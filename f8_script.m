%% Plots figure 8 for the manuscript.

% Contour plots for the full model with thinning and thickening.



%% 1. setup axes
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

clear
figure(8)
clf
tiledlayout(1,2,'TileSpacing','compact','padding','none')
set(gcf,'pos',[  57   436   736   324])
ax1 = nexttile;
ax2 = nexttile;


%% 2. define the vectors of beta, r_s and nu
% Load results from figure 4, just to get beta and r_s.
load savedResults/f4_results beta r_s_dim

% accumulation multiplier 
nu = [0.9 1.1];

% only compute up to value 15 for each parameter because the simulation
% fails for higher ones (just because the domain thins to zero or negative
% thickness).
betaMaxIndex = 15;
rsMaxIndex   = 15;
%% 3. run the full model
load rerun rerun
% rerun = 1;
if rerun 

    for kk = 1:length(nu)
        for ii  = 1:rsMaxIndex
            for jj = 1:betaMaxIndex
                p = FirnSetup5('scaleDuration',1,'simDuration',1,'beta',beta(jj),'r_s_dim',r_s_dim(ii),'nu',nu(kk));
                results_1(jj,ii,kk)  = FCM9b(p);
                zeta_830_full(jj,ii,kk) = results_1(jj,ii,kk).zeta830final;
                inner_loop_iter = betaMaxIndex-jj
            end
            outer_loop_iter = rsMaxIndex-ii
        end 
    end
    reduceResultsStruct
    save savedResults/f8_results p results_1 zeta_830_full r_s_dim beta betaMaxIndex rsMaxIndex
else
    load savedResults/f8_results
end


r_s = r_s_dim/p.r2_0;


%% 4. plot the results for thickening. 
contourf(ax1,beta(1:betaMaxIndex),r_s(1:rsMaxIndex),zeta_830_full(:,:,2),0:0.1:1,'ShowText','on')

%% 5. plot the results for thinning. 
contourf(ax2,beta(1:betaMaxIndex),r_s(1:rsMaxIndex),zeta_830_full(:,:,1),0:0.1:1,'ShowText','on')

%% 6. finish plots
% add panel labels
text(ax1,-0.15,1.04,'a','units','normalized','FontSize',20)
text(ax2,-0.15,1.04,'b','units','normalized','FontSize',20)

% make the color scale the same across panels and figures
caxis(ax1,[0 1])
caxis(ax2,[0 1])

% add axis labels
ylabel(ax1,'surface grain size, $r^2_s$')
xlabel(ax1,'baseline accumulation rate, $\beta$')
ylabel(ax2,'surface grain size, $r^2_s$')
xlabel(ax2,'baseline accumulation rate, $\beta$')

% add titles
title(ax1,'firn thickness, $z_{830}$; $\dot{h}>0$','FontSize',18)
title(ax2,'firn thickness, $z_{830}$; $\dot{h}<0$','FontSize',18)

% change size of tick labels
ax1 = axis_font_sizes(ax1,13,17);
ax2 = axis_font_sizes(ax2,13,17);


%% 6. print figure

print('-dpng','figures/f8.png')
