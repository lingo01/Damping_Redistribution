clear all; clc;
tic;

%% Initial Settings

% network settings
global nNet netY netB netG;
settingNetwork('Network_netY_IEEE39Hete');

% simulation settings
settingSimulation();
global LoadEquivTime LoadEquivTimeTmp LoadEquivDroop;
LoadEquivTime = 0.5;    %0.5
LoadEquivDroop = 50;    %50
LoadEquivTimeTmp = LoadEquivTime;

% element setting
global type CDInv Gen Load;
settingElement('ElementList_IEEE39Hete');

netY_orig = netY;   netG_orig = netG;   netB_orig = netB;
nNet_orig = nNet;

CDInv_orig = CDInv;
Load_orig = Load;
Gen_orig = Gen;
type_orig = type;

dynaDevPos = [04, 08, 12, 15, 20, 24, 30, 31, 32, 34, 35, 36, 37, 38, 39];

parallelAdmittance = 1e2j;

%% Dimension scaling

dynaDevNum = linspace(1, 1000/length(dynaDevPos), 25)' * ones( 1, length(dynaDevPos) );
dynaDevNum = round(dynaDevNum);

%%
for rr = 1 : size(dynaDevNum, 1)

    % Initialization 
    CDInv = CDInv_orig;
    Gen   = Gen_orig;
    netY = netY_orig;
    nNet = nNet_orig;
    type = type_orig;

    % admittance matrix reshape
    nNet_tmp = nNet + sum(dynaDevNum(rr,:));
    netY_tmp = zeros(nNet_tmp, nNet_tmp);
    netY_tmp(1:nNet, 1:nNet) = netY;

    for ii = 1:length(dynaDevPos)
        
        if ii == 1
            bias = nNet;
        else
            bias = nNet + sum(dynaDevNum(rr, 1:(ii-1)));
        end

        for jj = 1:dynaDevNum(rr, ii)
            netY_tmp(bias + jj, bias +jj) = -parallelAdmittance / dynaDevNum(rr,ii);
            netY_tmp(bias + jj, dynaDevPos(ii)) = parallelAdmittance / dynaDevNum(rr,ii);
            netY_tmp(dynaDevPos(ii), bias + jj) = parallelAdmittance / dynaDevNum(rr,ii);
            netY_tmp(dynaDevPos(ii), dynaDevPos(ii)) = netY_tmp(dynaDevPos(ii), dynaDevPos(ii)) - parallelAdmittance/ dynaDevNum(rr,ii);
        end   
    end

    % device control reshape
    for ii = 1:length(dynaDevPos)

        if ii == 1
            bias = nNet;
        else
            bias = nNet + sum(dynaDevNum(rr, 1:(ii-1)));
        end

        for jj = 1:dynaDevNum(rr, ii)
            type{bias + jj} = type{dynaDevPos(ii)};
            if strcmp(type{dynaDevPos(ii)}, "Gen")
                Gen{dynaDevPos(ii)} = Gen_orig{dynaDevPos(ii)};
                Gen{dynaDevPos(ii)}.D  = Gen{dynaDevPos(ii)}.D    ;%/ dynaDevNum(rr,ii);
                Gen{dynaDevPos(ii)}.xd = Gen{dynaDevPos(ii)}.xd   ;%* dynaDevNum(rr,ii);
                Gen{dynaDevPos(ii)}.xdp = Gen{dynaDevPos(ii)}.xdp ;%* dynaDevNum(rr,ii);
                Gen{dynaDevPos(ii)}.Pm = Gen{dynaDevPos(ii)}.Pm   / dynaDevNum(rr,ii);
                Gen{bias + jj} = Gen{dynaDevPos(ii)};
            elseif strcmp(type{dynaDevPos(ii)}, "CDInv")
                CDInv{dynaDevPos(ii)} = CDInv_orig{dynaDevPos(ii)};
                CDInv{dynaDevPos(ii)}.t1 = CDInv{dynaDevPos(ii)}.t1 ;%* dynaDevNum(rr,ii);
                CDInv{dynaDevPos(ii)}.D1 = CDInv{dynaDevPos(ii)}.D1 ;%* dynaDevNum(rr,ii);
                CDInv{dynaDevPos(ii)}.D2 = CDInv{dynaDevPos(ii)}.D2 ;%* dynaDevNum(rr,ii);
                CDInv{dynaDevPos(ii)}.Ps = CDInv{dynaDevPos(ii)}.Ps / dynaDevNum(rr,ii);
                CDInv{dynaDevPos(ii)}.Qs = CDInv{dynaDevPos(ii)}.Qs / dynaDevNum(rr,ii);
                CDInv{bias + jj} = CDInv{dynaDevPos(ii)};
            end
        end   

    end

    for ii = 1:length(dynaDevPos)
        type{dynaDevPos(ii)} = "Cnct";
    end

    % substitude
    netY = netY_tmp;
    netG = real(netY_tmp);
    netB = imag(netY_tmp);
    nNet = nNet_tmp;
    

    % steady state
    [state0, state0_calcflag] = funcSteady();

    mat = funMatrixReduce(state0, 1);

    eig_matT = eig(mat.T);
    eig_matT = real(eig_matT);
    eig_matT = sort(eig_matT, 'descend');
    eig_matT = eig_matT( 1 : (length(eig_matT)-2*39) );  

    devNum(rr) = sum(dynaDevNum(rr,:)) + (15 - length(dynaDevPos));

    totalDamping(rr)  = -sum(eig_matT);
    averageDamping(rr) = totalDamping(rr) / devNum(rr);

    fprintf('round %d: dev num = %d, total damping = %.4e, average damping = %.4e\n', rr, devNum(rr), totalDamping(rr), averageDamping(rr));
end

%%
figure(144);
subplot(3,2,1); hold on; grid on; box on; 
plot(devNum, totalDamping, '-o', 'LineWidth', 1.2);
xlabel('n_{total}', 'FontName', 'Times New Roman');
ylabel('Total damping', 'FontName', 'Times New Roman')
legend('non-uniform', 'uniform', 'FontName', 'Times New Roman')
title('(a)', 'FontName', 'Times New Roman')
xlim([0, 1000])

subplot(3,2,2); hold on; grid on; box on;
plot(devNum, totalDamping' ./ (dynaDevNum * [2;2;2;2;2;3;3;3;2;3;3;3;3;3;3]), '-o', 'LineWidth', 1.2);
xlabel('n_{total}', 'FontName', 'Times New Roman');
ylabel('Average damping', 'FontName', 'Times New Roman')
legend('non-uniform', 'uniform', 'FontName', 'Times New Roman')
title('(b)', 'FontName', 'Times New Roman')
xlim([0, 1000]);     ylim([0, Inf])

%%
toc;