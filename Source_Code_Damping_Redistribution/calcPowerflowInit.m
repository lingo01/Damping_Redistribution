% compute the initial values for the power flow solution
% Parameters:
%   x = [delta_1 V_1 delta_2 V_2 ...]

function xInit = calcPowerflowInit()
    global type Gen CDInv QDInv;
    global nNet;
	xInit = zeros(2*nNet,1);
    
    % solve the equilibrium point at each bus
    for ii = 1:nNet
        if strcmp(type{ii}, 'Gen')
            xInit(2*ii-1) = 0;
            xInit(2*ii)   = Gen{ii}.Ef;
        elseif strcmp(type{ii}, 'CDInv')
            xInit(2*ii-1) = CDInv{ii}.As;
            xInit(2*ii)   = CDInv{ii}.Vs;
        elseif strcmp(type{ii}, 'QDInv')
            xInit(2*ii-1) = QDInv{ii}.As;
            xInit(2*ii)   = QDInv{ii}.Vs;
        elseif strcmp(type{ii}, 'Load')
            xInit(2*ii-1) = 0;
            xInit(2*ii)   = 1.05;
        elseif strcmp(type{ii}, 'Cnct')
            xInit(2*ii-1) = 0;
            xInit(2*ii)   = 1.05;
        end
    end
end