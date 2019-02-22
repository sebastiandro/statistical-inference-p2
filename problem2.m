data_raw = readtable('./whales.txt', 'ReadVariableNames', false);
data = data_raw.Var1;
n = length(data); %% Number of counties

%% Task a'
bin_width = 0.28;
figure(1);
clf;
h = histogram(data, 20, 'BinWidth', bin_width);
hold on;
xlabel('t');
ylabel('Number of whales');

% Task b

% Get the first two moment estimates
mu1 = mean(data);
mu2 = mean(data.^2);

% Create alhpa and lambda estimates from moments
a_hat_mme = mu1^2/(mu2-mu1^2);
l_hat_mme = a_hat_mme/mu1;

gam = gampdf(0:0.01:max(data), a_hat_mme,l_hat_mme);
plot(0:0.01:max(data), gam*n*bin_width, 'LineWidth',1);

%% Task c
% Mostlikelihood gives us a_hat_mle = 1.5955
% The a_hat_mle can be found with function gamfit, but I found it using
% mathematica
a_hat_mle = 1.5954;
l_hat_mle = mean(data)/a_hat_mle;

%% Task d
gam_mle = gampdf(0:0.01:max(data), a_hat_mle,0.37);
plot(0:0.01:max(data), gam_mle*n*bin_width, 'LineWidth',1);
legend('Real Data', sprintf('Gam(%0.2f, %0.2f) MME', a_hat_mme,l_hat_mme), sprintf('Gam(%0.2f, %0.2f) MLE', a_hat_mle,l_hat_mle))

%% Task e
% Parametric bootstrap
N = 1500;
n = 200;

% Draw N x n samples
mme_gamma_dist_draws = random('Gamma', a_hat_mme, l_hat_mme, [N n]);
% From each sample of size n, take the first moment. 
mu1_mme_bootstrap = mean(mme_gamma_dist_draws, 2);
% From each sample of size n, take the second moment. 
mu2_mme_bootstrap = mean(mme_gamma_dist_draws .^2, 2);

% create Nx1 vector of alpha mme estimates
alphas = (mu1_mme_bootstrap.^2)./(mu2_mme_bootstrap - mu1_mme_bootstrap.^2);

% create Nx1 vector of lambda mme estimates
lambdas = alphas ./ mu1_mme_bootstrap;

figure(2);
clf;
histogram(alphas, 20)
title('MME Alpha Sampling Distribution')
xlabel('\alpha')
ylabel('Number of estimates')

figure(3);
clf;
histogram(lambdas, 20)
xlabel('\lambda')
ylabel('Number of estimates')
title('MME Lambda Sampling Distribution')

fprintf('Standard error for MME alpha sampling distribution: %0.4f\n', std(alphas) / sqrt(n));
fprintf('Standard error for MME lambda sampling distribution: %0.4f\n', std(lambdas) / sqrt(n))

%% Task f
N = 1500;
n = 200;

mle_lambas_alphas = zeros(N,2);
mle_gamma_dist_draws = random('Gamma', a_hat_mle, l_hat_mle, [N n]);
for i = 1:N
    sample = mle_gamma_dist_draws(i,:);
    % Here we use gamfit to find the MLE
    mle_lambas_alphas(i,:) = gamfit(sample);
end
% Plot alpha sampling distribution
figure(4);
histogram(mle_lambas_alphas(:,1),20)
title('MLE Alpha Sampling Distribution')
xlabel('\alpha')
ylabel('Number of estimates')

figure(5);
% Plot lambda sampling distribution
histogram(mle_lambas_alphas(:,2),20)
title('MLE Lambda Sampling Distribution')
xlabel('\lambda')
ylabel('Number of estimates')

mle_bootstrap_alpha_std = std(mle_lambas_alphas(:,1)) / sqrt(n);
mle_bootstrap_lambda_std = std(mle_lambas_alphas(:,2) / sqrt(n));

fprintf('Standard error for MLE alpha sampling distribution: %0.4f\n', mle_bootstrap_alpha_std);
fprintf('Standard error for MLE lambda sampling distribution: %0.4f\n', mle_bootstrap_lambda_std)

%% g
mle_params_sorted = sort(mle_lambas_alphas);

alpha_delta_lower = mle_params_sorted(round(N*0.05),1) - a_hat_mle;
alpha_delta_upper = mle_params_sorted(round(N*0.95),1) - a_hat_mle;

lambda_delta_lower = mle_params_sorted(round(N*0.05),2) - l_hat_mle;
lambda_delta_upper = mle_params_sorted(round(N*0.95),2) - l_hat_mle;

mle_alpha_interval = [a_hat_mle - alpha_delta_upper, a_hat_mle - alpha_delta_lower];
mle_lambda_interval = [l_hat_mle - lambda_delta_upper, l_hat_mle - lambda_delta_lower];