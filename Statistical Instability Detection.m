%% =========================================
% Phase 3
% Statistical Instability Detection using CUSUM
% =========================================

clear
clc
close all

%% Load NASA IMS Dataset

folder_path = "C:\Users\USER\Downloads\Projects\IMS\1st_test";

files = dir(folder_path);

count = 0;

for i = 1:length(files)

    name = files(i).name;

    if ~(files(i).isdir || endsWith(name,".m"))
        count = count + 1;
    end

end

signal = zeros(count,1);

index = 1;

for i = 1:length(files)

    name = files(i).name;

    if files(i).isdir || endsWith(name,".m")
        continue
    end

    file_path = fullfile(folder_path,name);

    data = load(file_path);

    vibration = data(:,1);

    signal(index) = rms(vibration);

    index = index + 1;

end

fprintf("Total measurements loaded: %d\n",length(signal));


%% Log Degradation State

g = log(signal);

%% Degradation Velocity

r = diff(g);

%% Adaptive Velocity Model

window = 50;

r_hat = movmean(r,window);

%% Residual Dynamics

epsilon = r - r_hat;

epsilon_norm = (epsilon - mean(epsilon)) / std(epsilon);

%% Sequential Instability Detection (CUSUM)

baseline = 200;

mu0 = mean(epsilon_norm(1:baseline));
sigma0 = std(epsilon_norm(1:baseline));

k = 0.5 * sigma0;
h = 5 * sigma0;

N = length(epsilon_norm);

S = zeros(N,1);

regime_points = zeros(N,1);
regime_count = 0;

S_current = 0;

for t = 2:N

    S_current = max(0,S_current + epsilon_norm(t) - mu0 - k);

    S(t) = S_current;

    if S_current > h

        regime_count = regime_count + 1;
        regime_points(regime_count) = t;

        S_current = 0;

    end

end

regime_points = regime_points(1:regime_count);

fprintf("\nDetected Regime Boundaries:\n")
disp(regime_points)

%% Plot CUSUM Statistic

figure
plot(S,'LineWidth',1.5)
hold on
yline(h,'r--','Threshold')

title("CUSUM Instability Detection")
xlabel("Time Index")
ylabel("CUSUM Statistic")

grid on
