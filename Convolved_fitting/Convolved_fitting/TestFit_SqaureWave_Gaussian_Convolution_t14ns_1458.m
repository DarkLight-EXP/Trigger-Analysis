% Import the CSV Data
data = readtable('/Users/gabbygelinas/Desktop/Masters/Convolved_fitting/csv/t14ns_1458');
timeDiff = data.timeDiff;
x_root = data.x_root;

% Define the Gaussian and Square Wave functions with amplitude parameter
gaussian = @(x, mu, sigma, A) A * exp(-0.5 * ((x - mu) / sigma).^2);
square_wave = @(x, width) double(abs(x) <= width / 2);

% Define the convolution model with amplitude parameter
conv_model = @(params, x) conv(gaussian(x, params(1), params(2), params(4)), ...
                               square_wave(x, params(3)), 'same');

% Estimate initial guesses for the parameters
[~, maxIdx] = max(timeDiff);
mu_initial = x_root(maxIdx); % Gaussian mean around the peak
sigma_initial = std(x_root); % Standard deviation based on the data spread
width_initial = range(x_root) / 15; % Width based on a fraction of the data range
amplitude_initial = max(timeDiff) / 35; % Initial guess for the amplitude, /13 is a good starting point

initial_guess = [mu_initial, sigma_initial, width_initial, amplitude_initial]; %initial guess array 

% Compute the convolution using the initial guess
initial_model = interp1(x_root, conv_model(initial_guess, x_root), x_root, 'linear', 'extrap');

% Plot the original data and the initial guess model (uncomment for later
% checks)

% figure;
% plot(x_root, timeDiff, 'b', 'DisplayName', 'Original Data');
% hold on;
% plot(x_root, initial_model, 'g', 'DisplayName', 'Initial Guess Model');
% legend;
% xlabel('x\_root');
% ylabel('timeDiff');
% title('Original Data and Initial Guess Model');
% grid on;

% Set the x-axis limits
xlim([min(x_root), max(x_root)]);

hold off;

% Define the objective function for fitting with amplitude parameter
objective = @(params) sum((interp1(x_root, conv_model(params, x_root), x_root, 'linear', 'extrap') - timeDiff).^2);

% Perform the fit using fminsearch
optimal_params = fminsearch(objective, initial_guess);

% Generate the fitted model
fitted_model = interp1(x_root, conv_model(optimal_params, x_root), x_root, 'linear', 'extrap');


%-------------------------------------------------------------------------

% Display the optimal parameters
disp('------ NEW FILE ------------');
disp('Optimal Parameters:');
disp(['Mu (mean): ', num2str(optimal_params(1))]);
disp(['Sigma (standard deviation): ', num2str(optimal_params(2))]);
disp(['Width: ', num2str(optimal_params(3))]);
disp(['Amplitude: ', num2str(optimal_params(4))]);

% Calculate the Full Width at Half Maximum (FWHM)
sigma_fitted = optimal_params(2);
FWHM = 2 * sqrt(2 * log(2)) * sigma_fitted;

disp(['Full Width at Half Maximum (FWHM): ', num2str(FWHM)]);

% Calculate the Residual Sum of Squares (RSS)
RSS = sum((timeDiff - fitted_model).^2);

% Calculate the Total Sum of Squares (TSS)
TSS = sum((timeDiff - mean(timeDiff)).^2);

% Calculate the Coefficient of Determination (R^2)
R_squared = 1 - (RSS / TSS);

% Calculate the Root Mean Square Error (RMSE)
RMSE = sqrt(mean((timeDiff - fitted_model).^2));

%Chi-Squared Fit .... Need further work to evaluate 
% Calculate the chi-squared statistic
%chi_squared = sum(((timeDiff - fitted_model) ./error_bars ).^2);
%chi_squared_red = sum(((timeDiff - fitted_model) ./error_bars ).^2)/(length(timeDiff) - 4);

% Display the fit quality metrics
disp('Fit Quality Metrics:');
disp(['Residual Sum of Squares (RSS): ', num2str(RSS)]);
disp(['Coefficient of Determination (R^2): ', num2str(R_squared)]);
disp(['Root Mean Square Error (RMSE): ', num2str(RMSE)]);
%disp(['Chi squared: ', num2str(chi_squared)]);
%disp(['Reduced Chi squared: ', num2str(chi_squared_red)]);



%-------------------------------------------------------------------------

% Create a figure with two subplots
figure;

% Plot the original data and the fitted model in the first subplot
subplot(2, 1, 1);
plot(x_root, timeDiff, 'b.', 'DisplayName', 'Original Data'); % Change line plot to scatter plot
hold on;
plot(x_root, fitted_model, 'r', 'DisplayName', 'Fitted Model');
legend;
xlabel('x\_root');
ylabel('timeDiff');
title('Original Data and Fitted Model');
grid on;

% Set the x-axis limits, this sets your x range
xlim([-2, 2]);

% Calculate the square root of timeDiff for error bars, poisson error
error_bars = sqrt(timeDiff);

% Plot error bars
errorbar(x_root, timeDiff, error_bars, 'b.', 'CapSize', 10);

hold off;

% Calculate the residuals
residuals = (timeDiff - fitted_model) ./ fitted_model;

% Plot the residuals vs fitted values in the second subplot
subplot(2, 1, 2);
plot(x_root, residuals, 'bo');
xlabel('x root');
ylabel('Normalized Residuals');
title('Normalized Residuals vs Fitted Values');
xlim([-2, 2]);
ylim([-1, 1]);
grid on;
refline(0, 0); % Add a reference line at y = 0

hold off;


% Set the parameters for the Gaussian and Square Wave
mu = optimal_params(1);
sigma = optimal_params(2);
width = optimal_params(3);
A = optimal_params(4);




% Define a function to plot the convolution (this was done to see if
% convolution function works...can uncomment for any checks)

% function plot_convolution(gaussian, square_wave, mu, sigma, width, A, x_range)
%     % Create a dense range of x values for plotting
%     x = linspace(min(x_range), max(x_range), 100); %can change number of points in the last entry
% 
%     % Compute the Gaussian and Square Wave
%     gauss_values = gaussian(x, mu, sigma, A);
%     square_wave_values = square_wave(x, width);
% 
%     % Perform the convolution
%     conv_values = conv(gauss_values, square_wave_values, 'same');
% 
%     % Plot the Gaussian, Square Wave, and their Convolution
%     figure;
%     subplot(3, 1, 1);
%     plot(x, gauss_values, 'r', 'DisplayName', 'Gaussian');
%     title('Gaussian Function');
%     xlabel('x');
%     ylabel('Amplitude');
%     grid on;
% 
%     subplot(3, 1, 2);
%     plot(x, square_wave_values, 'b', 'DisplayName', 'Square Wave');
%     title('Square Wave Function');
%     xlabel('x');
%     ylabel('Amplitude');
%     grid on;
% 
%     subplot(3, 1, 3);
%     plot(x, conv_values, 'm', 'DisplayName', 'Convolution');
%     title('Convolution of Gaussian and Square Wave');
%     xlabel('x');
%     ylabel('Amplitude');
%     grid on;
% 
%     % Display all legends
%     legend('show');
% end

% Plot the convolution (uncomment if wa
%plot_convolution(gaussian, square_wave, mu, sigma, width, A, x_root);
