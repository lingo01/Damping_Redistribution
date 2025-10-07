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

DevTimeScaling = [1, 1000/15];  

for rr = 1:length(DevTimeScaling)

for ii = 1:nNet
    if strcmp(type{ii}, 'CDInv')
        CDInv{ii}.D1 = CDInv_orig{ii}.D1 / DevTimeScaling(rr);
        CDInv{ii}.D2 = CDInv_orig{ii}.D2 / DevTimeScaling(rr);
    elseif strcmp(type{ii}, 'Gen') && ii ~= refPos
        Gen{ii}.M = Gen_orig{ii}.M * DevTimeScaling(rr);
        Gen{ii}.D = Gen_orig{ii}.D * DevTimeScaling(rr);
        Gen{ii}.xd = Gen_orig{ii}.xd / DevTimeScaling(rr);
        Gen{ii}.xdp = Gen_orig{ii}.xdp / DevTimeScaling(rr);
    end
end

% Steady state
[state0, state0_calcflag] = funcSteady();

% Matrix
mat = funMatrixReduce(state0);

% eigen-space
eig_matT = eig(mat.T);
eig_matT_uncouple = eig(mat.T_uncouple);

if rr == 1
    figure(144);     subplot(3,2,5);
    hold on;    box on;     grid on;
    plot(eig_matT_uncouple+1e-16j, 'o', 'LineWidth', 1.2);
    plot(eig_matT+1e-16j, '*', 'LineWidth', 1.2);
    title('(d)','FontName','Times New Roman');
    legend('eig(T_n)', 'eig(T)' ,'FontSize',6,'FontName','Times New Roman', 'Location', 'northwest');
    xlim([-1 0]);   
    hold off;
else
    figure(144);     subplot(3,2,6);
    hold on;    box on;     grid on;
    plot(eig_matT_uncouple+1e-16j, 'o', 'LineWidth', 1.2);
    plot(eig_matT+1e-16j, '*', 'LineWidth', 1.2);
    title('(e)','FontName','Times New Roman');
    legend('eig(T_n)', 'eig(T)' ,'FontSize',6,'FontName','Times New Roman', 'Location', 'northwest');
    xlim([-1 0]);   hold off;
end

end

%%
toc;