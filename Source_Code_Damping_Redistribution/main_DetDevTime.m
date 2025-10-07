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
LoadEquivTime_orig = LoadEquivTime;


%% Load scaling

DevTimeScaling = linspace(0.05, 2.0, 300);   

for rr = 1:length(DevTimeScaling)

% device time constant scaling
for ii = 1:nNet
    if strcmp(type{ii}, 'CDInv')
        CDInv{ii}.t1 = CDInv_orig{ii}.t1 * DevTimeScaling(rr);
    elseif strcmp(type{ii}, 'Gen') && ii ~= refPos
        Gen{ii}.M = Gen_orig{ii}.M * DevTimeScaling(rr);
    end
end

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


det_matT(rr) = det(mat.T);
det_matT_uncouple(rr) = det(mat.T_uncouple);
ratio_det(rr) = abs(prod(eig(mat.T)./eig(mat.T_uncouple)));

mod_index(rr) = norm(mat.t(2:end,2:end) / [mat.M, zeros(size(mat.M,1),size(mat.k,2)); zeros(size(mat.k,1),size(mat.M,2)), mat.k], inf)  * norm(mat.r+mat.g, inf) * norm(mat.D, 2) / min(eig(mat.A));
mod_index(rr) = 1/mod_index(rr);

ratio_fro_matT(rr) = norm(mat.T_couple,'fro')/norm(mat.T_uncouple,'fro');

fprintf('finished run rr = %d, scaling = %.16f, det ratio = %.6f\n', rr, DevTimeScaling(rr), ratio_det(rr))
end

%%
figure(200);
subplot(4,2,4);     box on; hold on; grid on;
yyaxis left;
plot(DevTimeScaling, ratio_det, 'LineWidth', 1.2);  ylim([0,1])
title('(d)','FontSize',10,'FontName','Times New Roman');
xlabel('\gamma_{\tau}','FontSize',10,'FontName','Times New Roman');
xlim([0.0, 2.0])

yyaxis right;
plot(DevTimeScaling, ratio_fro_matT, 'LineWidth', 1.2);
ylabel('||T_c||_F/||T_n||_F','FontSize',10,'FontName','Times New Roman');
ylim([0, 0.025]);

%%


%%
toc;