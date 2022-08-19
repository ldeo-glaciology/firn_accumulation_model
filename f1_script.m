%% %% Plots figure 1 for the manuscript and computes values 
%% for Table of scales and non-dimensional parameters

% The script computes scales and parameters by running the model function (porosity_model_v1.m), while
% passing it various values of the accumulation rate scale b_0 and the surface
% temperature T_s. It also passes porosity_model_v1.m 'justScalesParameterCalcs' = true, 
% so that the it knows not to actually run the model, just to compute the
% scales and output them. This is a slow way of doing this, but it ensures
% that the scales and parameter computed for the table and figure are the
% same as those used in the actual model.

% This version of the code uses the setup scripts from Rob Skarbek's
% method-of-lines code to get the values of the parameters and scales. 


% Modified in response to reviewer comments:
%    -- increased font size on figure axis labels and tick labels

%% 1. setup axes

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

clear 
figure(1)
clf
tiledlayout(1,2,'TileSpacing','compact','padding','normal')
set(gcf,'pos',[ 57   429   794   331])
ax1 = nexttile;
ax2 = nexttile;

%% 2. compute scales and non-dim parameters for Table 3 of the paper

% Define three values of b_0 and and three values of T_s
b_mpy = [1 0.1 0.01];
T_s_dim = 273.15 - [0 20 40];

% loop over these values
clear p temp
for rr = 1:3   
    p(rr) = FirnSetup5('b0_mpy',b_mpy(rr),'T_s_dim',T_s_dim(rr));     
end  


%% 3. collate the results into a table for writing to a latex table
parameterscale = {'$b_0$';'$T_s$';'$r_0$';'$z_0$';...
    '$w_0$';'$t_0$';'$\sigma_0$';'$A_r$';'$\delta$'};

description = {'accumulation scale [m a$^{-1}$]';'surface temperature [K]';'grain radius scale [m$^2$]'; 'vertical scale [m]';...
    'velocity scale [m a$^{-1}$]';'time scale [a]';'overburden stress scale [Pa]';...
    'Arthern number [-]';'grain size saturation ratio [-]'};

ii=1;
high = [p(ii).b0_mpy; p(ii).T_s; p(ii).r2_0; p(ii).z_0; ...
    p(ii).w_0*(365*24*60*60) ; p(ii).t_0/(365*24*60*60); p(ii).sigma_0; p(ii).ArthenNumber; p(ii).delta];

ii=2;
intermediate = [p(ii).b0_mpy; p(ii).T_s; p(ii).r2_0; p(ii).z_0; ...
    p(ii).w_0*(365*24*60*60) ; p(ii).t_0/(365*24*60*60); p(ii).sigma_0; p(ii).ArthenNumber; p(ii).delta];

ii=3;
low = [p(ii).b0_mpy; p(ii).T_s; p(ii).r2_0; p(ii).z_0; ...
    p(ii).w_0*(365*24*60*60) ; p(ii).t_0/(365*24*60*60); p(ii).sigma_0;p(ii).ArthenNumber; p(ii).delta];

% write a latex table
T = table(parameterscale, description, high, intermediate, low)
table2latex(T, 'table4.tex')


%% 4. compute values for the figure by looping over b_0 and T_0 simultaneously 
clear temp

b_mpy = 10.^[-2.:0.04:0];
Tm = 273.15;
T_s_dim = Tm - [0:1:40];
load rerun rerun
if rerun
    for rr = 1:length(b_mpy)
        for ss = 1:length(T_s_dim)
            p = FirnSetup5('b0_mpy',b_mpy(rr),'T_s_dim',T_s_dim(ss));
            delta(rr,ss) = p.delta;
            Ar(rr,ss) = p.ArthenNumber;
            Pe(rr,ss) = p.PecletNumber;
        end
        length(b_mpy)-rr
    end
    save savedResults/f1_results p delta Ar Pe
else
    load savedResults/f1_results p delta Ar Pe
end

%% 5. plot the figure
%%% panel A: alpha
contourf(ax1,T_s_dim,b_mpy,(Ar),'ShowText','on')
set(ax1,'YDir','normal','YScale','log')
ylabel(ax1,'surface accmulation scale, $b_0$ [m yr$^{-1}]$','interpreter','latex');
xlabel(ax1,'surface temperature, $T_s$ [K]','interpreter','latex');
title_a = title(ax1,'compaction number, $\alpha$','FontSize',15);
title_b.Position(2) = 1.1;
ylim(ax1,[1e-2, 1])
xlim(ax1,[233 273.15])

%%% panel B: delta
contourf(ax2,T_s_dim,b_mpy,log10(delta),'ShowText','on')
set(ax2,'YDir','normal','YScale','log')
xlabel('surface temperature, $T_s$ [K]','interpreter','latex');
title_b = title(ax2,'logarithm of grain size saturation ratio, log$_{10}\delta$','FontSize',15);
title_b.Position(2) = 1.1;
ylim([1e-2, 1])
xlim([233 273.15])

% increase size of tick labels
ax1.XAxis.FontSize = 12;
ax1.XLabel.FontSize = 15;
ax1.YAxis.FontSize = 12;
ax1.YLabel.FontSize = 15;
ax2.XAxis.FontSize = 12;
ax2.XLabel.FontSize = 15;
ax2.YAxis.FontSize = 12;
ax2.YLabel.FontSize = 15;

% add panel labels
text(ax1,-0.18,1.06,'(a)','units','normalized','FontSize',20)
text(ax2,-0.18,1.06,'(b)','units','normalized','FontSize',20)

%% 5. print figure
print('-dpng','figures/f1.png') 
