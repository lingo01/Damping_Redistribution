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

DauxCtrlScaling = 2e-3;
auxCtrlScaling = -0.5:DauxCtrlScaling:1.5;   

for rr = 1:length(auxCtrlScaling)

% Steady state
[state0, state0_calcflag] = funcSteady();

% Matrix
mat = funMatrixReduce(state0);

mat.aux = zeros(size(mat.T));
numV1 = size(mat.M,1);
numV2 = size(mat.k,1);
numV = numV1 + numV2;

tmp1 = [inv(mat.M), zeros(numV1, numV2); zeros(numV2, numV1), inv(mat.k)];
tmp2 = inv(mat.t);
for ii = 1:size(mat.k)
    if tmp1(ii,ii) == 1/ (LoadEquivTime/LoadEquivDroop)
        tmp1(ii,ii) = 0;
    end
end
for ii = 1:size(mat.t,1)
    if tmp2(ii,ii) == 1/LoadEquivTime
        tmp2(ii,ii) = 0;
    end
end

mat.aux( (numV1+1):(2*numV1+numV2), (2*numV1+numV2+2):end ) = - 1 * tmp1 * auxCtrlScaling(rr) * eye(numV);

mat.aux( (2*numV1+numV2+2):(3*numV1+numV2+1), 1:numV1 ) =  0.01 * auxCtrlScaling(rr) * eye(numV1);
mat.aux( (2*numV1+numV2+2):(3*numV1+numV2+1), (numV1+1):2*numV1 ) =  5 * auxCtrlScaling(rr) * eye(numV1);
mat.aux( (3*numV1+numV2+2):end, (2*numV1+1):(2*numV1+numV2) ) = 1 * auxCtrlScaling(rr) * eye(numV2);

mat.aux( (2*numV1+numV2+1):end, 1:(2*numV1+numV2) ) = - tmp2 * mat.aux( (2*numV1+numV2+1):end, 1:(2*numV1+numV2) );

mat.Taux = mat.T + mat.aux;

tmp = sort(real(eig(mat.Taux)));    tmp = tmp(end:-1:1);
maxeig_matTaux(rr,1) = tmp(1);
tmp1 = tmp(find(tmp<tmp(1)));
maxeig_matTaux(rr,2) = tmp1(1);
tmp2 = tmp(find(tmp<tmp1(1)));
maxeig_matTaux(rr,3) = tmp2(1);


if rr == 486% 401 % 701 % 486
    eig2 = eig(mat.Taux);
elseif rr == round((0.2+0.5)/DauxCtrlScaling) + 1
    eig3 = eig(mat.Taux);
elseif auxCtrlScaling(rr) == 0
    eig1 = eig(mat.Taux);
end

ratio_det(rr) = abs(prod(eig(mat.Taux)./eig(mat.T)));

fprintf('finished run rr = %d, scaling = %.3f, det ratio = %.6f\n', rr, auxCtrlScaling(rr), ratio_det(rr))
end

%% plots
close all;

figure(134); 

subplot(3,2,[1,2]);  
hold on;  box on; grid on;
h1 = plot(auxCtrlScaling, ratio_det, 'LineWidth', 1.2);
h3 = fill(auxCtrlScaling(1)-DauxCtrlScaling+DauxCtrlScaling*[min(find(ratio_det>=1)), max(find(ratio_det>=1)), max(find(ratio_det>=1)), min(find(ratio_det>=1))], [0.85,0.85, 1.05, 1.05], 'c'); 
set(h3, 'FaceColor', [0 0.4470 0.7410]);
set(h3, 'FaceAlpha', 0.2);  set(h3, 'EdgeAlpha', 0)
xlabel('\sigma', 'FontName', 'Times New Roman')
ylabel('||T''''|/|T||', 'FontName', 'Times New Roman')
xlim([min(auxCtrlScaling), max(auxCtrlScaling)]);
ylim([0.85, 1.05])
title('(a)', 'FontName', 'Times New Roman')


subplot(3,2,[3,4]); 
hold on; box on; grid on;
h2 = plot(auxCtrlScaling, maxeig_matTaux(:,1), 'LineWidth', 1.2, "Color", [0.85,0.33,0.10]);
h4 = plot(auxCtrlScaling, maxeig_matTaux(:,2), '--', 'LineWidth', 1.2, "Color", [0.85,0.33,0.10]);
h5 = plot(auxCtrlScaling, maxeig_matTaux(:,3), ':', 'LineWidth', 1.2, "Color", [0.85,0.33,0.10]);
h3 = fill(auxCtrlScaling(1)-DauxCtrlScaling+DauxCtrlScaling*[min(find(ratio_det>=1)), max(find(ratio_det>=1)), max(find(ratio_det>=1)), min(find(ratio_det>=1))], [-0.5,-0.5, 0,0], 'c'); 
set(h3, 'FaceColor', [0 0.4470 0.7410]);
set(h3, 'FaceAlpha', 0.2);  set(h3, 'EdgeAlpha', 0)

xlabel('\sigma', 'FontName', 'Times New Roman')
ylabel('Re(\lambda)', 'FontName', 'Times New Roman')
xlim([min(auxCtrlScaling), max(auxCtrlScaling)])
ylim([-0.5, 0])
legend([h2,h4,h5], 'max Re(\lambda)', 'max Re_{(2)}(\lambda)', 'max Re_{(3)}(\lambda)', 'FontName', 'Times New Roman', 'Location', 'northeast')
title('(b)', 'FontName', 'Times New Roman')

%
a1 = sort(real(eig1)); a1 = a1(end:-1:1);  a1 = a1(find(a1>-0.5));
a3 = sort(real(eig3)); a3 = a3(end:-1:1);  a3 = a3(find(a3>-0.5));

subplot(3,2,5);
hold on;  grid on; box on;
plot(eig1, 's', 'LineWidth', 1.2);
plot(eig3, '+', 'LineWidth', 1.2);
legend('\sigma=0', '\sigma=0.2', 'FontName', 'Times New Roman', 'FontSize', 8, 'Location', 'northwest');  
xlim([-0.5 -0]);
title('(c)', 'FontName', 'Times New Roman')
xlabel('Re(\lambda)', 'FontName', 'Times New Roman')
ylabel('Im(\lambda)', 'FontName', 'Times New Roman')
for ii = 1:length(a1)
    tmp = find(real(eig1)==a1(ii));
    if imag(eig1(tmp(1))) ~= 0 
        text(a1(ii)+0.01, imag( eig1( tmp(1) ) ), ['\textcircled{', num2str(ii), '}'], 'FontName', 'Times New Roman', 'FontSize', 8, 'Interpreter','latex');
    else
        text(a1(ii)-0.01, imag( eig1( tmp(1) ) )-2.5, ['\textcircled{', num2str(ii), '}'], 'FontName', 'Times New Roman', 'FontSize', 8, 'Interpreter','latex');
    end
end

subplot(3,2,6);
tmp = min(length(a1), length(a3));
tmp = a3(1:tmp)  - a1(1:tmp);  
bar(tmp)
title('(d)', 'FontName', 'Times New Roman')
xlabel('eigenvalue ID', 'FontName', 'Times New Roman')
ylabel('Re(\lambda'''')-Re(\lambda)', 'FontName', 'Times New Roman');



%%
toc;