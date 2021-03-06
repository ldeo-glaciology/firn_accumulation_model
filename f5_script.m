%% Plots figure 5 for the manuscript.

% Plot the sensitivity of firn thickness to accumulation rate and how this 
% changes with surface grain size. It does this by loading the results of
% running the model many times with many combinations of beta and r_s
% (using F4_beta_evolving_r_justFullModel_v3.m), then using a linear fit to
% get the the rate of change of z_830 with beta.


% Modified in response to reviewer comments:
%    -- increased font size on figure axis labels and tick labels

%% 1. setup axes
clear 
figure(5)
clf

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

tiledlayout(1,1,'TileSpacing','none','padding','none')
set(gcf,'pos',[57   501   332   259])
ax1 = nexttile;

%% 2. load the parameter-search results
% load the results of running the full model many times for many values
% of r_s and beta - to rerun these results run
% F4_beta_evolving_r_justFullModel_v3.m with rerun = 1
load savedResults/f4_results


%% 3. compute dzeta830db with a linear least-squares fit.   
for ii =1:size(results_1,2)
    r_s(ii) = results_1(1,ii).p.r2_s_dim/results_1(1,ii).p.r2_0;
    zeta_830(:,ii) = [results_1(:,ii).zeta830final];
    temp = polyfit(beta',zeta_830(:,ii),1);
    dzeta830db(ii,1) = temp(1);
end

%% 4. plot dzeta830db as a function of r_s
plot(ax1,r_s, dzeta830db,'-r','LineWidth',1)    % automatically stops at N =18 because at higher r_s the zeta_820 arrays contain NaNs
hold on
scatter(ax1,r_s,dzeta830db,20,'r','filled')

%% 5. finish figure
xlabel(ax1,'surface grain size, $r^2_s$')
ylabel(ax1,{'mean dependence of firn', 'thickness on $\beta$'})

% change size of tick labels
ax1 = axis_font_sizes(ax1,13,17);



%% 6. print figure
print('-dpng','figures/f5.png')


