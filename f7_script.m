%% Plots figure 7 for the manuscript.

% Produces three contour plot of firn thickness, zeta_830, as a function of
% accumulation rate, beta, and surface grain size, r_s. Each contour plot corresponds
% to a different value of the stress exponent, n (2, 3 and 4). 

% Modified in response to reviewer comments:
%    -- increased font size on figure axis labels and tick labels

%% 1. setup axes
set(groot, 'defaulttextinterpreter','latex');
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
% Load results from figure 4, just to get beta and r_s.
load savedResults/f4_results

% the stress exponent
n = [1 2 3 4];

%% 3. run the full model for every combination of beta, r_s and n
load rerun rerun
% rerun = 1;
if rerun 
for kk = 1:length(n)
    for ii  = 1:length(r_s_dim)
        for jj = 1:length(beta)
            p = FirnSetup5('beta',beta(jj),'r_s_dim',r_s_dim(ii),'n',n(kk));
            results_1(jj,ii,kk)  = FCM9b(p);
            zeta_830_full(jj,ii,kk) = results_1(jj,ii,kk).zeta830final;           
        inner_loop_iter = length(beta)-jj
        end
        middle_loop_iter = length(r_s_dim)-ii
    end
    outer_loop_iter = length(n)-kk
end

% save the results for use later
reduceResultsStruct % remove all but the final timestep in each of the variables, using this short script. 
save savedResults/f7_results results_1 zeta_830_full r_s beta
end

%% 4. reload the results of this loop, to save time when replotting
if ~rerun
    load savedResults/f7_results
end

%% 5. plot three contour maps, n=2, n=3 and n=4
contourf(ax1,beta,r_s,zeta_830_full(:,:,2)','ShowText','on')
contourf(ax2,beta,r_s,zeta_830_full(:,:,3)','ShowText','on')
contourf(ax3,beta,r_s,zeta_830_full(:,:,4)','ShowText','on')

%% 6. finish figure
% add panel labels
text(ax1,-0.15,0.98,'a','units','normalized','FontSize',20)
text(ax2,-0.15,0.98,'b','units','normalized','FontSize',20)
text(ax3,-0.15,0.98,'c','units','normalized','FontSize',20)

% add acis labels
ylabel(ax1,'surface grain size, $r^2_s$')
xlabel(ax1,'accumulation rate, $\beta$')
ylabel(ax2,'surface grain size, $r^2_s$')
xlabel(ax2,'accumulation rate, $\beta$')
ylabel(ax3,'surface grain size, $r^2_s$')
xlabel(ax3,'accumulation rate, $\beta$')

% make the color scale the same across panels and figures
caxis(ax1,[0 1])
caxis(ax2,[0 1])
caxis(ax3,[0 1])

% add titles
title(ax1,'$n$ = 2','FontSize',18)
title(ax2,'$n$ = 3','FontSize',18)
title(ax3,'$n$ = 4','FontSize',18)

% change size of tick labels
ax1 = axis_font_sizes(ax1,13,17);
ax2 = axis_font_sizes(ax2,13,17);
ax3 = axis_font_sizes(ax3,13,17);


%% 7. print figure
print('-dpng','figures/f7.png')
