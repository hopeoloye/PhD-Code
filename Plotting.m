% Assuming Areas is defined as:
Areas = {'L_A1','R_A1','L_A4','R_A4','L_44','R_44'}; % Example areas

% Generate parameter names for field 'A'
parameter_names = {};
for i = 1:length(Areas)
    for j = 1:length(Areas)
        parameter_names{end+1} = ['A(', Areas{i}, ' -> ', Areas{j}, ')'];
    end
end

% Define your true values (if available, otherwise use zeros)
num_params = length(PEB.Ep); % Number of estimated parameters
true_values = zeros(1, num_params); % Example true values, replace if known

% Plot EP with error bars and true values for the loaded PEB results
plot_EP_with_true_values(PEB.Ep, PEB.Cp, true_values, parameter_names);

% Define your true values (example values, replace with actual true values)
true_values = [0.5, -0.3, 0.2, 0.4, -0.1, 0.3, -0.4, 0.1, 0.2, 0.5]; % Example true values
parameter_names = {'Param1', 'Param2', 'Param3', 'Param4', 'Param5', 'Param6', 'Param7', 'Param8', 'Param9', 'Param10'}; % Replace with actual parameter names

% Plot EP with error bars and true values for the loaded PEB results
plot_EP_with_true_values(PEB.Ep, PEB.Cp, true_values, parameter_names);

% Additional PEB results plotting if needed
% load(['PEB_HC_vs_B2_no_covars_v' Ver '.mat']);
% plot_EP_with_true_values(PEB.Ep, PEB.Cp, true_values, parameter_names);

% load(['PEB_HC_vs_B3_no_covars_v' Ver '.mat']);
% plot_EP_with_true_values(PEB.Ep, PEB.Cp, true_values, parameter_names);

% Create a heatmap of the design matrix M.X
figure;
heatmap(M.X, 'GridVisible', 'off');

% Add title and labels
title('Design Matrix M.X');
xlabel('Covariates');
ylabel('Subjects');


% Create an image plot of the design matrix M.X
figure;
imagesc(M.X);

% Add colorbar to indicate the scale
colorbar;

% Add title and labels
title('Design Matrix M.X');
xlabel('Covariates');
ylabel('Subjects');

% Set x-axis and y-axis ticks
xticks(1:size(M.X, 2));
yticks(1:size(M.X, 1));

% Label the ticks if you have specific names for covariates and subjects
% xticklabels({'Intercept', 'Covariate 1', 'Covariate 2', '...'});
% yticklabels({'Subject 1', 'Subject 2', 'Subject 3', '...'});



% Assuming M is already defined and contains the design matrix

figure;
imagesc(M.X);
colorbar;
xlabel('Covariates and Group');
ylabel('Subjects');
title('Design Matrix (M.X)'); figure

% conflict
