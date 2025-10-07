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

DevTimeScaling = [0.5, 1.0, 1.5, 2.5];  
titleSerial = {'(a)', '(b)', '(c)', '(d)'};

for rr = 1:length(DevTimeScaling)
    
for ii = 1:nNet
    if strcmp(type{ii}, 'CDInv')
        CDInv{ii}.D2 = CDInv_orig{ii}.D2 * DevTimeScaling(rr);
    elseif strcmp(type{ii}, 'Gen') && ii ~= refPos
        Gen{ii}.xd = Gen_orig{ii}.xd * DevTimeScaling(rr);
        Gen{ii}.xdp = Gen_orig{ii}.xdp * DevTimeScaling(rr);
    end
end

% Steady state
[state0, state0_calcflag] = funcSteady();

% Matrix
mat = funMatrixReduce(state0);

% eigen-space
eig_matT = eig(mat.T);
eig_matT_uncouple = eig(mat.T_uncouple);

% figure(30);     subplot(2,2,rr);
figure(198);     subplot(6,2,rr+4);
hold on;    box on;     grid on;
plot(eig_matT_uncouple, 'o', 'LineWidth', 1.2);
plot(eig_matT, '*', 'LineWidth', 1.2);
title(titleSerial{rr},'FontSize',10,'FontName','Times New Roman');
xlim([-1 0]);   hold off;

fprintf('det(T)/det(T_uncouple) = %.4f\n', abs(prod(eig(mat.T)./eig(mat.T_uncouple))));

end

%%
toc;