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

loadScaling = (1-0.75*logspace(0,-5,300)) * 1.28972/0.6;     % 1.28972/0.6 

for rr = 1:length(loadScaling)

% device time constant scaling
for ii = 1:nNet
    if strcmp(type{ii}, 'Load')
        Load{ii}.Ps = Load_orig{ii}.Ps * loadScaling(rr);
        Load{ii}.Qs = Load_orig{ii}.Qs * loadScaling(rr);
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

ratio_det(rr) = abs(prod(eig_matT(:,rr)./eig(mat.T_uncouple)));

mod_index(rr) = norm(mat.t(2:end,2:end) / [mat.M, zeros(size(mat.M,1),size(mat.k,2)); zeros(size(mat.k,1),size(mat.M,2)), mat.k], inf) * norm(mat.r+mat.g, inf) * norm(mat.D, 2) / min(eig(mat.A));
mod_index(rr) = 1/mod_index(rr);

ratio_fro_matT(rr) = norm(mat.T_couple,'fro')/norm(mat.T_uncouple,'fro');

fprintf('finished run rr = %d, scaling = %.16f, det ratio = %.6f\n', rr, loadScaling(rr), ratio_det(rr))
end

%%
figure(100);    hold on;
subplot(3,1,1); plot(loadScaling, maxeig_matT, 'LineWidth', 1.2);   title('max eig(T)')
subplot(3,1,2); plot(loadScaling, maxeig_matTa, 'LineWidth', 1.2);  title('max eig(Ta)')
subplot(3,1,3); plot(loadScaling, maxeig_matTv, 'LineWidth', 1.2);  title('max eig(Tv)')
hold off;

figure(200);
subplot(4,2,1);     box on; hold on; grid on;
yyaxis left;
plot(loadScaling, ratio_det, 'LineWidth', 1.2);   ylim([0,1]);
title('(a)','FontSize',10,'FontName','Times New Roman');
xlabel('\rho_{load}','FontSize',10,'FontName','Times New Roman');   
ylabel('|det(T)/det(T_n)|','FontSize',10,'FontName','Times New Roman')
xlim([0.5,2.2])
xline(loadScaling(end), '--', "Color", "#A2142F", 'LineWidth', 1.2);              
text(loadScaling(end)-1.5,0.9,'voltage collapse \leftarrow', 'FontName', 'Times New Roman')

yyaxis right;
plot(loadScaling, ratio_fro_matT, 'LineWidth', 1.2);
ylim([0, 0.025]);

%%
toc;