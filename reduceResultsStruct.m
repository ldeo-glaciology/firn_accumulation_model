

%% reduce size of results structure
% by removing everything except the last column of the large x-t arrays

for jj = 1:length(results_1(:))
    results_1(jj).Rho = results_1(jj).Rho(:,end);
    results_1(jj).W = results_1(jj).W(:,end);
    results_1(jj).Temp = results_1(jj).Temp(:,end);
    results_1(jj).GrainSize = results_1(jj).GrainSize(:,end);
    results_1(jj).Age = results_1(jj).Age(:,end);
    results_1(jj).Phi = results_1(jj).Phi(:,end);
    results_1(jj).Sigma = results_1(jj).Sigma(:,end);
    results_1(jj).Depth = results_1(jj).Depth(:,end);
end