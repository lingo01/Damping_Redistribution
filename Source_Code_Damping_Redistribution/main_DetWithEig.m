clear all; clc;
tic;

%% Initial Settings

% network settings
global nNet netB;
settingNetwork('Network_netY_IEEE39Hete');

% simulation settings
global refPos LoadEquivTime LoadEquivDroop;
settingSimulation();
global LoadEquivTimeTmp
LoadEquivTimeTmp = LoadEquivTime;

% element setting
global type CDInv Gen Load;
settingElement('ElementList_IEEE39Hete');

CDInv_orig = CDInv;
Load_orig = Load;
Gen_orig = Gen;

%% Load scaling

sampleNum = 10000;
loadScaling =  rand(sampleNum, 1) * (2.15-0.5) + 0.5;
DevCtrlScaling =  rand(sampleNum, 1) * (2.5-0.05) + 0.05;
DevTimeScaling =  rand(sampleNum, 1) * (1-0.05) + 0.05;

ratio_det = zeros(sampleNum, 1) * NaN;
max_eig_Tn = zeros(sampleNum, 1) * NaN;
max_eig_T  = zeros(sampleNum, 1) * NaN;

for rr = 1:length(loadScaling)

% device time constant scaling
for ii = 1:nNet
    if strcmp(type{ii}, 'Load')
        Load{ii}.Ps = Load_orig{ii}.Ps * loadScaling(rr);
        Load{ii}.Qs = Load_orig{ii}.Qs * loadScaling(rr);
    elseif strcmp(type{ii}, 'CDInv')
        CDInv{ii}.t1 = CDInv_orig{ii}.t1 * DevTimeScaling(rr);
        CDInv{ii}.D2 = CDInv_orig{ii}.D2 * DevCtrlScaling(rr);
    elseif strcmp(type{ii}, 'Gen') && ii ~= refPos
        Gen{ii}.M = Gen_orig{ii}.M * DevTimeScaling(rr);
        Gen{ii}.xd = Gen_orig{ii}.xd * DevCtrlScaling(rr);
        Gen{ii}.xdp = Gen_orig{ii}.xdp * DevCtrlScaling(rr);
    end
end

% Steady state
[state0, state0_calcflag] = funcSteady();

if ~state0_calcflag
    % Matrix
    mat = funMatrixReduce(state0);
    
    % eigen-space
    max_eig_Tn(rr) = max(real(eig(mat.T_uncouple)));
    max_eig_T(rr) = max(real(eig(mat.T)));
    ratio_det(rr) = abs(det(mat.T)/det(mat.T_uncouple));

    slope_actual(rr) = (ratio_det(rr) - 1) / (max_eig_T(rr) - max_eig_Tn(rr));
    slope_approx(rr) = 1 / max_eig_Tn(rr);
end

fprintf('finished run rr = %d\n', rr)
end


max_eigen_shift = max_eig_T-max_eig_Tn;
max_eigen_shift_ratio_det = ratio_det;


%% plot

figure(200);
subplot(4,2,7); 
% subplot(2,3,[1,2]);  
box on; hold on; grid on;
scatter(max_eig_Tn, ratio_det, 05, 'LineWidth', 0.3);
scatter(max_eig_T,  ratio_det, 05, 'LineWidth', 0.3);       
% xlim([-0.2, 0])
legend('max R(\lambda^{(n)})','max R(\lambda)','FontSize',9,'FontName','Times New Roman')
title('(e)','FontSize',10,'FontName','Times New Roman');
xlabel('max R(\lambda) or max R(\lambda^{(n)})','FontSize',10,'FontName','Times New Roman');   
ylabel('|det(T)/det(T_n)|','FontSize',10,'FontName','Times New Roman')
ylim([1e-3, 1]);

subplot(4,2,8); 
box on; hold on; grid on;
scatter(max_eig_T-max_eig_Tn, ratio_det, 05, 'LineWidth', 0.3);
title('(f)','FontSize',10,'FontName','Times New Roman');
xlabel('max R(\lambda) - max R(\lambda^{(n)})','FontSize',10,'FontName','Times New Roman');   
ylim([1e-3, 1]);

approx_line_x = 0:0.001:0.1; approx_line_y = -1/0.1234 * (approx_line_x) + 0.82;
plot(approx_line_x, approx_line_y, '--', 'LineWidth', 0.5)

legend('actual shift','approx shift','FontSize',9,'FontName','Times New Roman')

%%
toc;