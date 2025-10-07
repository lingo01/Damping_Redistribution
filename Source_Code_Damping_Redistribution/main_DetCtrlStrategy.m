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
global type CDInv QDInv Gen Load;
settingElement('ElementList_IEEE39Hete');

CDInv_orig = CDInv;
QDInv_orig = QDInv;
Load_orig = Load;
Gen_orig = Gen;
type_orig = type;

[state0, state0_calcflag] = funcSteady();
power0 = calcPower(state0);
P0 = real(power0);    Q0 = imag(power0);

mat = funMatrixReduce(state0);
ratio_det0 = abs( prod (eig(mat.T) ./ eig(mat.T_uncouple) ) );

% sum_cap = sum(P0(find(P0>=0)));
sum_cap = sum(abs(P0));
GFL_cap = 0;
for ii = 1:nNet
    if strcmp(type{ii},"Load")
        % if Load{ii}.Ps >= 0
            GFL_cap = GFL_cap + abs(Load{ii}.Ps);
        % end
    end
end
ratio_cap0 = GFL_cap / sum_cap;

%% sampling numbers of GFL replacement trajectory
nSamp = 100;

%% Load scaling

unsortSerial = [4, 8, 12, 15, 20, 24, 30, 31, 32, 34, 35, 36, 37, 38]; % position of static loads in IEEE 39-bus system

sortSerial = zeros(nSamp, length(unsortSerial));
for ii = 1:nSamp
    idx = randperm(length(unsortSerial));
    sortSerial(ii, :) = unsortSerial(idx);
end

ratio_det = NaN * zeros(size(sortSerial));
ratio_cap = NaN * zeros(size(sortSerial));

for rr = 1:size(sortSerial,1)

    Gen = Gen_orig;
    CDInv = CDInv_orig;
    QDInv = QDInv_orig;
    Load = Load_orig;

    type = type_orig;

    for kk = 1:size(sortSerial, 2)

        Load{sortSerial(rr,kk)}.Ps = P0(sortSerial(rr,kk));
        Load{sortSerial(rr,kk)}.Qs = Q0(sortSerial(rr,kk));
        type{sortSerial(rr,kk)} = "Load";

        % Steady state
        [state0, state0_calcflag] = funcSteady();

        if state0_calcflag == 1
            break;
        end

        % state dynamic matrix
        mat = funMatrixReduce(state0);
        eig_matT = eig(mat.T);
        eig_matTn = eig(mat.T_uncouple);

        if max(real(eig_matT)) > 0
            break;
        end

        ratio_det(rr,kk) = abs( prod (eig_matT ./ eig_matTn ) );

        if kk >= 2
            if ratio_det(rr,kk) > ratio_det(rr,kk-1)
                ratio_det(rr,kk) = NaN;
                break;
            end
        end

        GFL_cap = 0;
        for ii = 1:nNet
            if strcmp(type{ii},"Load")
                GFL_cap = GFL_cap + abs(Load{ii}.Ps);
            end
        end
        ratio_cap(rr,kk) = GFL_cap / sum_cap;

        fprintf('finished run rr = %d, kk = %d, det = %.6f\n', rr, kk, ratio_det(rr,kk))
    end
end

ratio_det = [ones(size(sortSerial,1), 1)*ratio_det0, ratio_det];
ratio_cap = [ones(size(sortSerial,1), 1)*ratio_cap0, ratio_cap];

%%
figure(200);    
subplot(4,2,[5,6]);
box on; grid on; hold on;
for ii = 1:size(ratio_det,1)
    plot(ratio_cap(ii,:), ratio_det(ii,:), '-*', 'LineWidth', 0.5);
end
xlabel('GFL penetration', 'FontSize', 10, 'FontName', 'Times New Roman');
ylabel('|det(T)/det(Tn)|', 'FontSize', 10, 'FontName', 'Times New Roman')
yline(0.4, '--', 'Color', 'r', 'LineWidth', 1.5);
title('(g)','FontSize',10,'FontName','Times New Roman');
xlim([0.25,0.95]);
hold off;


%%
toc;