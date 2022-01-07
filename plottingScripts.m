%% Scripts that plot the figures for the manuscript

%%% 
% To re-run all the simulations set rerun = 1
% if rerun = 0,  in some cases saved results will be reloaded to save
% time. 
clear all
rerun = 0
save rerun rerun
close all


%% figure 1 and the scale and parameters table:

f1_script 
 
%% figure 2 (comparison between full model and ODE model results)

f2_script    
    
%% figure 3 (dependence on accumulation rate with no grain size evolution for the full model and ODE model)

f3_script    
      
%% figure 4 (dependence on accumulation rate with grain size evolution for the full model)

f4_script    
    
%% figure 5 (dependence on accumulation rate with grain size evolution for the full model 2: gradient plot)

f5_script    
     
%% figure 6    (contour plots of firn thickness for full model and two ODE models)

f6_script    

%% figure 7 (contour plots for three nonlinear rheologys, n ~= 1)
  
f7_script    

%% figure 8 (contour plots of thickening and thinning experiments). 

f8_script    

