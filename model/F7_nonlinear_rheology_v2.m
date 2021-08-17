%% Plots figure 7 for the manuscript.

% Produces three contour plot of firn thickness, zeta_830, as a function of
% accumulation rate, beta, and surface grain size, r_s. Each contour plot corresponds
% to a different value of the stress exponent, n (2, 3 and 4). 

% This version of the code uses the setup scripts from Rob Skarbek's
% method-of-lines code.

%% 1. setup axes
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

clear
figure(7)
clf
tiledlayout(1,3,'TileSpacing','none','padding','none')
set(gcf,'pos',[ 57         436        1069         324])
ax1 = nexttile;
ax2 = nexttile;
ax3 = nexttile;

%% 2. define the vectors of beta, r_s and n
% load results from figure 4 to get beta and r_s
load results_for_Figure4_4

% the stress exponent
n = [1 2 3 4];

%% 3. run the full model for every combination of beta, r_s and n
load rerun rerun
if rerun 
for kk = 1:length(n)
    for ii  = 1:length(r_s_dim)
        for jj = 1:length(beta)
            p = FirnSetup3('beta',beta(jj),'r_s_dim',r_s_dim(ii),'n',n(kk));
            results_1(jj,ii,kk)  = FCM9(p);
            zeta_830_full(jj,ii,kk) = results_1(jj,ii,kk).zeta830final;           
        end
        ii
    end
end

% save the results for use later
reduceResultsStruct % remove all but the final timestep in each of the variables, using this short script. 
save results_for_nonlinearrheology_plots_v2 results_1 zeta_830_full r_s beta
end

%% 4. reload the results of this loop, to save time when replotting
% (to rerun the results, uncomment the cell above and comment out the load
% command below, then rerun the whole script.) 
if ~rerun
    load results_for_nonlinearrheology_plots_v2
end

%% 5. plot three contour maps, n=2, n=3 and n=4
contourf(ax1,beta,r_s,zeta_830_full(:,:,2)','ShowText','on')
contourf(ax2,beta,r_s,zeta_830_full(:,:,3)','ShowText','on')
contourf(ax3,beta,r_s,zeta_830_full(:,:,4)','ShowText','on')


%% 6. finish figure
text(ax1,-0.11,0.98,'a','units','normalized','FontSize',20)
text(ax2,-0.11,0.98,'b','units','normalized','FontSize',20)
text(ax3,-0.11,0.98,'c','units','normalized','FontSize',20)

ylabel(ax1,'surface grain size, $r^2_s$','FontSize',15)
xlabel(ax1,'accumulation rate, $\beta$','FontSize',15)
ylabel(ax2,'surface grain size, $r^2_s$','FontSize',15)
xlabel(ax2,'accumulation rate, $\beta$','FontSize',15)
ylabel(ax3,'surface grain size, $r^2_s$','FontSize',15)
xlabel(ax3,'accumulation rate, $\beta$','FontSize',15)

caxis(ax1,[0 1])
caxis(ax2,[0 1])
caxis(ax3,[0 1])

title(ax1,'$n$ = 2','FontSize',15)
title(ax2,'$n$ = 3','FontSize',15)
title(ax3,'$n$ = 4','FontSize',15)


%% 7. print figure
print('-dpng','F7_nonlinear_contours_v2.png')
