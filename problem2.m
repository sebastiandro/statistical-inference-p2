data_raw = readtable('../data/whales.txt', 'ReadVariableNames', false);
data = data_raw.Var1;
N = length(data); %% Number of counties

%% Task a
clf;
histogram(data, 20);
hold on

% Task b

% Get the first two moment estimates
mu1 = mean(data);
mu2 = mean(data.^2);

% Create alhpa and lambda estimates from moments
l_hat = mu1/(mu2-mu1^2);    
a_hat = l_hat*mu1;

gam = gampdf(0:0.01:5, a_hat,l_hat);
plot(0:0.01:5, gam*60);
legend('Real Data', sprintf('Gam(%0.2f, %0.2f)', a_hat,l_hat))