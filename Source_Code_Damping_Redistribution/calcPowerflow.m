% compute the power mismatch at each bus
% Parameters:
%   x = [delta_1 V_1 delta_2 V_2 ...]

function f = calcPowerflow(x)   
    global type Gen CDInv QDInv Load;
    global nNet;

	f = zeros(2*nNet,1);
    
    % solve the steady-state power at each bus
    Psteady = zeros(nNet,1);
    Qsteady = zeros(nNet,1);
    for ii = 1:nNet
        if strcmp(type{ii}, 'Gen')
            Psteady(ii) = Gen{ii}.Pm;
            Qsteady(ii) = (Gen{ii}.Ef-x(ii*2))*x(ii*2)/(Gen{ii}.xd-Gen{ii}.xdp);
        elseif strcmp(type{ii}, 'CDInv')
            Psteady(ii) = CDInv{ii}.Ps - 1/CDInv{ii}.D1*(x(ii*2-1)-CDInv{ii}.As);
            Qsteady(ii) = CDInv{ii}.Qs - 1/CDInv{ii}.D2*(x(ii*2)  -CDInv{ii}.Vs);
        elseif strcmp(type{ii}, 'QDInv')
            Psteady(ii) = QDInv{ii}.Ps - 1/QDInv{ii}.D1*(x(ii*2-1)-QDInv{ii}.As);
            us = QDInv{ii}.Vs + QDInv{ii}.D2*QDInv{ii}.Qs/QDInv{ii}.Vs;
            Qsteady(ii) = -1/QDInv{ii}.D2 * x(ii*2) * (x(ii*2)-us);
        elseif strcmp(type{ii}, 'Load')
            Psteady(ii) = Load{ii}.Ps;
            Qsteady(ii) = Load{ii}.Qs;
        elseif strcmp(type{ii}, 'Cnct')
            Psteady(ii) = 0;
            Qsteady(ii) = 0;
        end
    end
    
    S = calcPower(x);   Pout = real(S); Qout = imag(S);
    for ii = 1:nNet
        f(2*ii-1) = Pout(ii) - Psteady(ii);
        f(2*ii)   = Qout(ii) - Qsteady(ii);
    end
end