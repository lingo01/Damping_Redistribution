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

dynaDevNumSum = 500 * rand(4000, 1);

dynaDevNum = zeros(length(dynaDevNumSum), length(dynaDevPos));
for ii = 1 : length(dynaDevNumSum)
    tmp = rand(1, length(dynaDevPos));
    dynaDevNum(ii, :) = dynaDevNumSum(ii, 1) * tmp / sum(tmp);
end
dynaDevNum = round(dynaDevNum) + 1;

% dynaDevNum = 100/length(dynaDevPos) * ones(1, 1) * ones( 1, length(dynaDevPos) );
% dynaDevNum = round(dynaDevNum);

% dynaDevNum = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] + 2;


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

    % figure(431); plot(state0(1:2:end), '-o', 'LineWidth', 1.2)
    % figure(432); plot(state0(2:2:end), '-o', 'LineWidth', 1.2)

    mat = funMatrixReduce(state0, 1);

    eig_matT = eig(mat.T);
    eig_matT = real(eig_matT);
    eig_matT = sort(eig_matT, 'descend');
    eig_matT = eig_matT( 1 : (length(eig_matT)-2*39) );  

    % diag_matT = diag(mat.T); 
    % diag_matT = diag_matT(find(diag_matT>-1e3));

    devNum(rr) = sum(dynaDevNum(rr,:)) + (15 - length(dynaDevPos));
    devNum_indiv(rr,:) = dynaDevNum(rr,:);
    % totalDamping(rr)  = -sum(diag_matT);
    totalDamping(rr)  = -sum(eig_matT);
    averageDamping(rr) = totalDamping(rr) / devNum(rr);

    fprintf('round %d: dev num = %d, total damping = %.4e, average damping = %.4e\n', rr, devNum(rr), totalDamping(rr), averageDamping(rr));
end

%%
figure(144);
subplot(3,2,1); hold on; grid on; box on; 
scatter(devNum, totalDamping, 10, 'o', 'LineWidth', 0.3);
xlabel('n_{total}');
ylabel('Total Damping')

subplot(3,2,2); hold on; grid on; box on;
scatter(devNum, totalDamping'./(devNum_indiv*[2;2;2;2;2;3;3;3;2;3;3;3;3;3;3]), 10, 'o', 'LineWidth', 0.3);
xlabel('n_{total}');
ylabel('Average Damping')


%%
toc;