clear all; clc;
tic;

%% Initial Settings

% network settings
global nNet netB;
settingNetwork('Network_netY_IEEE39Hete');

% simulation settings
global refPos LoadEquivTime LoadEquivDroop;
settingSimulation();

LoadEquivTime  = 0.5;
global LoadEquivTimeTmp
LoadEquivTimeTmp = LoadEquivTime;

% element setting
global type CDInv Gen Load;
settingElement('ElementList_IEEE39Hete');

CDInv_orig = CDInv;
Load_orig = Load;
Gen_orig = Gen;

netB_orig = netB;

%% Load scaling

netScaling = linspace(0.5, 3, 300);

for rr = 1:length(netScaling)

% device time constant scaling
netB = netB_orig * netScaling(rr);

% Steady state
[state0, state0_calcflag] = funcSteady();

% Matrix
mat = funMatrixReduce(state0);

% eigen-space
eig_matT(:,rr) = eig(mat.T);
eig_matTa(:,rr) = eig(mat.Ta);
eig_matTv(:,rr) = eig(mat.Tv);

maxeig_matT(rr) = max(real(eig_matT(:,rr)));
maxeig_matTa(rr) = max(real(eig_matTa(:,rr)));
maxeig_matTv(rr) = max(real(eig_matTv(:,rr)));

maxeig_matCut(rr) = max(real(eig( mat.Tv-mat.Tva*inv(mat.Ta)*mat.Tav )));

maxeig_matA(rr) = max(real(eig(mat.A)));
maxeig_matC1(rr) = max(real(eig(mat.C1)));
maxeig_matC2(rr) = max(real(eig(mat.C2)));

ratio_det(rr) = abs(prod(eig_matT(:,rr)./eig(mat.T_uncouple)));

mod_index(rr) = norm(mat.t(2:end,2:end) / [mat.M, zeros(size(mat.M,1),size(mat.k,2)); zeros(size(mat.k,1),size(mat.M,2)), mat.k], inf) * norm(mat.r+mat.g, inf) * norm(mat.D, 2) / min(eig(mat.A));
mod_index(rr) = 1/mod_index(rr);

ratio_fro_matT(rr) = norm(mat.T_couple,'fro')/norm(mat.T_uncouple,'fro');

fprintf('finished run rr = %d, scaling = %.16f, det ratio = %.6f, min damp = %.6f\n', rr, netScaling(rr), ratio_det(rr), -maxeig_matT(rr))
end

%% plot 

figure(200);
subplot(4,2,2);     box on; hold on; grid on;
yyaxis left;
plot(netScaling, ratio_det, 'LineWidth', 1.2);   ylim([0,1]);   xlim([0.5 4])
title('(b)','FontSize',10,'FontName','Times New Roman');
xlabel('\rho_{net}','FontSize',10,'FontName','Times New Roman');   

yyaxis right;
plot(netScaling, ratio_fro_matT, 'LineWidth', 1.2);
ylabel('||T_c||_F/||T_n||_F','FontSize',10,'FontName','Times New Roman');
ylim([0, 0.025]);

%%
toc;